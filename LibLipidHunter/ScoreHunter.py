# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2017  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
# LipidHunter is Dual-licensed
#     For academic and non-commercial use: `GPLv2 License` Please read more information by the following link:
#         [The GNU General Public License version 2] (https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
#     For commercial use:
#         please contact the SysMedOs_team by email.
# Please cite our publication in an appropriate form.
# Ni, Zhixu, Georgia Angelidou, Mike Lange, Ralf Hoffmann, and Maria Fedorova.
# "LipidHunter identifies phospholipids by high-throughput processing of LC-MS and shotgun lipidomics datasets."
# Analytical Chemistry (2017).
# DOI: 10.1021/acs.analchem.7b01126
#
# For more info please contact:
#     SysMedOs_team: oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
#     Developer Georgia Angelidou georgia.angelidou@uni-leipzig.de

from __future__ import division

import json
import math

import pandas as pd

from LibLipidHunter.ParallelFunc import ppm_calc_para, ppm_window_para, pr_window_calc_para
from LibLipidHunter.IsotopeHunter import IsotopeHunter
from LibLipidHunter.ScoreGenerator import ScoreGenerator
from LibLipidHunter.PanelPlotter import plot_spectra


def get_lipid_info(param_dct, fa_df, checked_info_df, checked_info_groups, core_list, usr_weight_df,
                   usr_key_frag_df, core_spec_dct, xic_dct):

    print('===> ==> -->Start to hunt for lipids!')

    usr_lipid_type = param_dct['lipid_type']
    charge_mode = param_dct['charge_mode']
    output_folder = param_dct['img_output_folder_str']
    usr_ms2_threshold = param_dct['ms2_th']
    usr_ms1_precision = param_dct['ms_ppm'] * 1e-6
    usr_ms2_ppm = param_dct['ms2_ppm']
    usr_ms2_info_th = param_dct['ms2_infopeak_threshold']
    usr_hg_info_th = param_dct['ms2_hginfopeak_threshold']

    usr_isotope_score_filter = param_dct['isotope_score_filter']
    usr_rank_score_filter = param_dct['rank_score_filter']

    img_typ = param_dct['img_type']
    img_dpi = param_dct['img_dpi']
    usr_fast_isotope = param_dct['fast_isotope']

    hunter_start_time_str = param_dct['hunter_start_time']
    isotope_hunter = IsotopeHunter()

    score_calc = ScoreGenerator(param_dct, usr_weight_df, usr_key_frag_df, usr_lipid_type,
                                checked_info_df, ion_charge=charge_mode, ms2_ppm=usr_ms2_ppm)

    tmp_df = pd.DataFrame()

    for group_key in core_list:
        _subgroup_df = checked_info_groups.get_group(group_key)
        usr_spec_info_dct = core_spec_dct[group_key]
        _samemz_se = _subgroup_df.iloc[0, :].squeeze()

        _usr_ms2_rt = _samemz_se['scan_time']
        _usr_ms2_dda_rank = _samemz_se['DDA_rank']
        _usr_ms2_scan_id = _samemz_se['scan_number']
        _usr_mz_lib = _samemz_se['Lib_mz']
        _usr_formula = _samemz_se['FORMULA']
        _usr_charge = _samemz_se['Ion']
        if _usr_charge == '[M-H]-':
            _usr_formula_charged = _samemz_se['[M-H]-_FORMULA']
            _usr_ms2_pr_mz = _samemz_se['[M-H]-_MZ']
        else:
            _usr_formula_charged = _samemz_se['[M-H]-_FORMULA']
            _usr_ms2_pr_mz = _samemz_se['Lib_mz']
        _ms1_pr_i = usr_spec_info_dct['ms1_i']
        _ms1_pr_mz = usr_spec_info_dct['ms1_mz']
        _ms1_df = usr_spec_info_dct['ms1_df']
        _ms2_df = usr_spec_info_dct['ms2_df']

        print(_usr_ms2_rt, _ms1_pr_mz, _usr_formula_charged)

        # use the max threshold from abs & relative intensity settings
        if 'i' in _ms2_df.columns.tolist():
            _ms2_max_i = _ms2_df['i'].max()
            ms2_threshold = max(usr_ms2_threshold, _ms2_max_i * usr_ms2_info_th)
            _score_ms2_df = _ms2_df.query('i > %f' % ms2_threshold)
        else:
            _ms2_df = pd.DataFrame()
            _score_ms2_df = pd.DataFrame()
        if _ms1_pr_mz > 0.0 and _ms1_df.shape[0] > 0 and _ms2_df.shape[0] > 0 and _ms1_pr_i > 0.0:

            print('>>> >>> >>> >>> Best PR on MS1: %f' % _ms1_pr_mz)

            isotope_score_info_dct = isotope_hunter.get_isotope_score(_ms1_pr_mz, _ms1_pr_i,
                                                                      _usr_formula_charged, _ms1_df,
                                                                      ms1_precision=usr_ms1_precision,
                                                                      isotope_number=2,
                                                                      only_c=usr_fast_isotope,
                                                                      score_filter=usr_isotope_score_filter)

            isotope_score = isotope_score_info_dct['isotope_score']

            print('isotope_score: %f' % isotope_score)
            if isotope_score >= usr_isotope_score_filter:
                print('>>> isotope_check PASSED! >>> >>> >>>')
                print('>>> >>> >>> >>> Entry Info >>> >>> >>> >>> ')
                _samemz_se.set_value('MS1_obs_mz', _ms1_pr_mz)
                _exact_ppm = 1e6 * (_ms1_pr_mz - _usr_mz_lib) / _usr_mz_lib
                _samemz_se.set_value('ppm', _exact_ppm)
                _samemz_se.set_value('abs_ppm', abs(_exact_ppm))

                for _i_abbr, _r_abbr in _subgroup_df.iterrows():
                    _usr_abbr_bulk = _r_abbr['BULK_ABBR']
                    print('Check_proposed_structure:', _usr_abbr_bulk)

                    matched_checker, obs_info_dct = score_calc.get_rankscore(fa_df, checked_info_df, _usr_abbr_bulk,
                                                                             charge_mode, _usr_mz_lib, _score_ms2_df,
                                                                             usr_lipid_type, ms2_ppm=usr_ms2_ppm,
                                                                             rankscore_filter=usr_rank_score_filter)

                    obs_info_df = obs_info_dct['INFO']
                    rank_score = obs_info_df['RANK_SCORE'].tolist()

                    if matched_checker > 0:
                        print('Rank_score: --> passed', rank_score)
                        print(obs_info_df[['BULK_ABBR', 'DISCRETE_ABBR', 'RANK_SCORE', 'scan_time']])

                        # format abbr. for file names
                        _save_abbr_bulk = _usr_abbr_bulk
                        _save_abbr_bulk = _save_abbr_bulk.replace(r'(', r'[')
                        _save_abbr_bulk = _save_abbr_bulk.replace(r')', r']')
                        _save_abbr_bulk = _save_abbr_bulk.replace(r'<', r'[')
                        _save_abbr_bulk = _save_abbr_bulk.replace(r'>', r']')
                        _save_abbr_bulk = _save_abbr_bulk.replace(r':', r'-')
                        _save_abbr_bulk = _save_abbr_bulk.replace(r'@', r'-')
                        _save_abbr_bulk = _save_abbr_bulk.replace('\\', r'_')
                        _save_abbr_bulk = _save_abbr_bulk.replace(r'/', r'_')

                        img_name_core = ('/%.4f_rt%.3f_DDAtop%.0f_scan%.0f_%s.%s'
                                         % (_usr_ms2_pr_mz, _usr_ms2_rt, _usr_ms2_dda_rank,
                                            _usr_ms2_scan_id, _save_abbr_bulk, img_typ)
                                         )
                        img_name = (output_folder +
                                    r'/LipidHunter_Results_Figures_%s'
                                    % hunter_start_time_str + img_name_core)

                        isotope_checker, isotope_score, img_n = plot_spectra(_usr_abbr_bulk,
                                                                             _samemz_se,
                                                                             xic_dct,
                                                                             obs_info_dct,
                                                                             usr_spec_info_dct,
                                                                             isotope_score_info_dct,
                                                                             _usr_formula_charged,
                                                                             _usr_charge,
                                                                             save_img_as=img_name,
                                                                             img_type=img_typ,
                                                                             dpi=img_dpi,
                                                                             ms1_precision=
                                                                             usr_ms1_precision
                                                                             )

                        print('==> check for output -->')

                        obs_info_df['Proposed_structures'] = _usr_abbr_bulk
                        obs_info_df['Bulk_identification'] = _usr_abbr_bulk
                        obs_info_df['Formula_neutral'] = _usr_formula
                        obs_info_df['Formula_ion'] = _usr_formula_charged
                        obs_info_df['Charge'] = _usr_charge
                        obs_info_df['MS1_obs_mz'] = _ms1_pr_mz
                        obs_info_df['MS1_obs_i'] = '%.2e' % float(_ms1_pr_i)
                        obs_info_df['Lib_mz'] = _usr_mz_lib
                        obs_info_df['MS2_scan_time'] = _usr_ms2_rt
                        obs_info_df['DDA#'] = _usr_ms2_dda_rank
                        obs_info_df['MS2_PR_mz'] = _usr_ms2_pr_mz
                        obs_info_df['Scan#'] = _usr_ms2_scan_id
                        # obs_info_df['#Specific_peaks'] = (target_frag_count + target_nl_count)
                        # obs_info_df['#Contaminated_peaks'] = (other_frag_count + other_nl_count)
                        obs_info_df['ppm'] = _exact_ppm

                        # if any IO error while writing img output
                        if img_n == '-2':
                            obs_info_df['img_name'] = '%s-2.%s' % (img_name_core[1:-4], img_typ)
                        else:
                            obs_info_df['img_name'] = img_name_core[1:]

                        tmp_df = tmp_df.append(obs_info_df)

    if tmp_df.shape[0] > 0:
        print('Size of the identified LPP_df %i, %i' % (tmp_df.shape[0], tmp_df.shape[1]))
        tmp_df.reset_index(drop=True, inplace=True)
        tmp_df.index += 1
    else:
        print('!! Size of the identified LPP_df == 0')
        tmp_df = 'error'

    return tmp_df
