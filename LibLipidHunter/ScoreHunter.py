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

from LibLipidHunter.ScoreGenerator import ScoreGenerator
from LibLipidHunter.PanelPlotter import plot_spectra
from LibLipidHunter.IsotopeHunter import IsotopeHunter


def get_lipid_info(param_dct, checked_info_df, checked_info_groups, core_list, usr_fa_def_df, usr_weight_df,
                   usr_key_frag_df, usr_scan_info_df, ms1_xic_mz_lst, core_spec_dct, xic_dct, target_ident_lst):
    tmp_idx = 1
    tmp_df = pd.DataFrame()

    usr_lipid_type = param_dct['lipid_type']
    charge_mode = param_dct['charge_mode']
    usr_vendor = param_dct['vendor']

    output_folder = param_dct['img_output_folder_str']

    usr_dda_top = param_dct['dda_top']
    usr_ms2_threshold = param_dct['ms2_th']
    usr_ms1_precision = param_dct['ms_ppm'] * 1e-6
    usr_ms2_precision = param_dct['ms2_ppm'] * 1e-6
    usr_hg_precision = param_dct['hg_ppm'] * 1e-6

    usr_overallscore_filter = param_dct['score_filter']
    usr_ms2_info_th = param_dct['ms2_infopeak_threshold']
    usr_hg_info_th = param_dct['ms2_hginfopeak_threshold']

    usr_fast_isotope = param_dct['fast_isotope']

    usr_sn_ratio = param_dct['sn_ratio']
    usr_isotope_score_filter = param_dct['isotope_score_filter']
    usr_rank_score_filter = param_dct['rank_score_filter']
    usr_msp_score_filter = param_dct['msp_score_filter']
    usr_fp_score_filter = param_dct['fp_score_filter']
    usr_snr_score_filter = param_dct['snr_score_filter']

    img_typ = param_dct['img_type']
    img_dpi = param_dct['img_dpi']
    usr_fast_isotope = param_dct['fast_isotope']

    # use the SNR equation SNR = 20 * log10(signal/noise)
    # snr_score = 20 * math.log10((signal_sum_i / noise_sum_i))
    # set s/n == 20 --> SNR_SCORE = 100
    # default 10.4795 = 100 / (20 * math.log10(3)) --> 10.4795
    if usr_sn_ratio == 3 or usr_sn_ratio == 0:
        usr_amp_factor = 10.4795
    else:
        usr_amp_factor = 100 / (20 * math.log10(usr_sn_ratio))

    if usr_fast_isotope is True:
        isotope_score_mode = '(Fast mode)'
    else:
        isotope_score_mode = ''

    hunter_start_time_str = param_dct['hunter_start_time']
    isotope_hunter = IsotopeHunter()

    score_calc = ScoreGenerator(param_dct, usr_fa_def_df, usr_weight_df, usr_key_frag_df, usr_lipid_type,
                                checked_info_df, ion_charge=charge_mode, ms2_ppm=param_dct['ms2_ppm'])

    for group_key in core_list:
        _subgroup_df = checked_info_groups.get_group(group_key)
        usr_spec_info_dct = core_spec_dct[group_key]
        _samemz_se = _subgroup_df.iloc[0, :].squeeze()
        _usr_ms2_pr_mz = _samemz_se['MS2_PR_mz']
        _usr_ms2_rt = _samemz_se['scan_time']
        _usr_ms2_dda_rank = _samemz_se['DDA_rank']
        _usr_ms2_scan_id = _samemz_se['scan_number']
        _usr_mz_lib = _samemz_se['Lib_mz']
        _usr_formula = _samemz_se['FORMULA_NEUTRAL']
        _usr_formula_charged = _samemz_se['Formula']
        _ms1_pr_i = usr_spec_info_dct['ms1_i']
        _ms1_pr_mz = usr_spec_info_dct['ms1_mz']
        _ms1_df = usr_spec_info_dct['ms1_df']
        _ms2_df = usr_spec_info_dct['ms2_df']

        _usr_charge = _samemz_se['Ion']
        print(_usr_ms2_rt, _ms1_pr_mz, _usr_formula_charged)

        # use the max threshold from abs & relative intensity settings
        if 'i' in _ms2_df.columns.tolist():
            _ms2_max_i = _ms2_df['i'].max()
            ms2_threshold = max(usr_ms2_threshold, _ms2_max_i * usr_ms2_info_th)
            _score_ms2_df = _ms2_df.query('i > %f' % ms2_threshold)
        else:
            _ms2_df = pd.DataFrame()
        if _ms1_pr_mz > 0.0 and _ms1_df.shape[0] > 0 and _ms2_df.shape[0] > 0 and _ms1_pr_i > 0.0:

            print('>>> >>> >>> >>> Best PR on MS1: %f' % _ms1_pr_mz)

            isotope_score_info_dct = isotope_hunter.get_isotope_score(_ms1_pr_mz, _ms1_pr_i,
                                                                      _usr_formula_charged, _ms1_df,
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

                other_signals_dct = score_calc.get_fa_possibilities(_usr_mz_lib, _usr_charge, _score_ms2_df)
                if len(other_signals_dct.keys()) > 0:
                    for _i_abbr, _r_abbr in _subgroup_df.iterrows():
                        _usr_abbr_bulk = _r_abbr['Abbreviation']
                        print('Check_proposed_structure:', _usr_abbr_bulk)

                        match_info_dct, matched_checker = score_calc.get_rankscore(_usr_abbr_bulk, charge_mode,
                                                                                   _usr_mz_lib, _score_ms2_df,
                                                                                   other_signals_dct, _ms2_max_i)
                        rank_score = match_info_dct['Rank_score']

                        if matched_checker > 0 and rank_score > usr_rank_score_filter:
                            print('Rank_score: %.f --> passed' % rank_score)

                            # snr_score, usr_sn_ratio, noise_df, snr_i_info = score_calc.get_snr_score(
                            #     match_info_dct,
                            #     specific_check_dct,
                            #     amplify_factor=
                            #     usr_amp_factor,
                            #     use_fp=
                            #     use_fp)
                            #
                            # if usr_sn_ratio >= 0 and snr_score >= usr_snr_score_filter:
                            #
                            #     overall_score = sum([rank_score, snr_score, isotope_score]) / 3
                            #     overall_score = round(overall_score, 1)
                            #     print('overall_score', overall_score)
                            #
                            #     if overall_score >= usr_overallscore_filter:
                            #
                            #         match_info_dct['Isotope_score'] = isotope_score
                            #         match_info_dct['SNR_score'] = snr_score
                            #         match_info_dct['Overall_score'] = overall_score
                            #
                            #         # format abbr. for file names
                            #         _save_abbr_bulk = _usr_abbr_bulk
                            #         _save_abbr_bulk = _save_abbr_bulk.replace(r'(', r'[')
                            #         _save_abbr_bulk = _save_abbr_bulk.replace(r')', r']')
                            #         _save_abbr_bulk = _save_abbr_bulk.replace(r'<', r'[')
                            #         _save_abbr_bulk = _save_abbr_bulk.replace(r'>', r']')
                            #         _save_abbr_bulk = _save_abbr_bulk.replace(r':', r'-')
                            #         _save_abbr_bulk = _save_abbr_bulk.replace(r'@', r'-')
                            #         _save_abbr_bulk = _save_abbr_bulk.replace('\\', r'_')
                            #         _save_abbr_bulk = _save_abbr_bulk.replace(r'/', r'_')
                            #
                            #         img_name_core = ('/%.4f_rt%.3f_DDAtop%.0f_scan%.0f_%s.%s'
                            #                          % (_usr_ms2_pr_mz, _usr_ms2_rt, _usr_ms2_dda_rank,
                            #                             _usr_ms2_scan_id, _save_abbr_bulk, img_typ)
                            #                          )
                            #         img_name = (output_folder +
                            #                     r'/LPPtiger_Results_Figures_%s'
                            #                     % hunter_start_time_str + img_name_core)
                            #
                            #         isotope_checker, isotope_score, img_n = plot_spectra(_usr_abbr_bulk,
                            #                                                              _samemz_se,
                            #                                                              xic_dct,
                            #                                                              match_info_dct,
                            #                                                              usr_spec_info_dct,
                            #                                                              specific_check_dct,
                            #                                                              isotope_score_info_dct,
                            #                                                              _usr_formula_charged,
                            #                                                              _usr_charge,
                            #                                                              save_img_as=img_name,
                            #                                                              img_type=img_typ,
                            #                                                              dpi=img_dpi,
                            #                                                              ms1_precision=
                            #                                                              usr_ms1_precision,
                            #                                                              msp_info=_msp_df,
                            #                                                              obs_fp=obs_fp_lst,
                            #                                                              missed_fp=missed_fp_lst,
                            #                                                              snr_i_info=snr_i_info
                            #                                                              )
                            #
                            #         print('==> check for output -->')
                            #
                            #         if _ms1_pr_i > 0 and isotope_checker == 0 \
                            #                 and isotope_score > usr_isotope_score_filter:
                            #
                            #             if 'OTHER_FRAG' in specific_check_dct.keys():
                            #                 other_frag_df = specific_check_dct['OTHER_FRAG']
                            #                 other_frag_count = other_frag_df.shape[0]
                            #             else:
                            #                 other_frag_count = 0
                            #             if 'OTHER_NL' in specific_check_dct.keys():
                            #                 other_nl_df = specific_check_dct['OTHER_NL']
                            #                 other_nl_count = other_nl_df.shape[0]
                            #             else:
                            #                 other_nl_count = 0
                            #             if 'TARGET_FRAG' in specific_check_dct.keys():
                            #                 target_frag_df = specific_check_dct['TARGET_FRAG']
                            #                 target_frag_count = target_frag_df.shape[0]
                            #                 target_frag_col_lst = target_frag_df.columns.tolist()
                            #                 for _frag_abbr in target_frag_col_lst:
                            #                     if _frag_abbr not in ['mz', 'i', 'LABEL', 'CLASS']:
                            #                         for _i, _f_se in target_frag_df.iterrows():
                            #                             if _f_se['LABEL'] == _frag_abbr:
                            #                                 match_info_dct[_frag_abbr] = _f_se[_frag_abbr]
                            #                                 if _frag_abbr not in target_ident_lst:
                            #                                     target_ident_lst.append(_frag_abbr)
                            #             else:
                            #                 target_frag_count = 0
                            #             if 'TARGET_NL' in specific_check_dct.keys():
                            #                 target_nl_df = specific_check_dct['TARGET_NL']
                            #                 target_nl_count = target_nl_df.shape[0]
                            #                 target_nl_col_lst = target_nl_df.columns.tolist()
                            #                 for _nl_abbr in target_nl_col_lst:
                            #                     if _nl_abbr not in ['mz', 'i', 'LABEL', 'CLASS']:
                            #                         for _i, _n_se in target_nl_df.iterrows():
                            #                             if _n_se['LABEL'] == _nl_abbr:
                            #                                 match_info_dct[_nl_abbr] = _n_se[_nl_abbr]
                            #                                 if _nl_abbr not in target_ident_lst:
                            #                                     target_ident_lst.append(_nl_abbr)
                            #             else:
                            #                 target_nl_count = 0
                            #
                            #             match_info_dct['Proposed_structures'] = _usr_abbr_bulk
                            #             match_info_dct['Bulk_identification'] = _usr_abbr_bulk
                            #             match_info_dct['Formula_neutral'] = _usr_formula
                            #             match_info_dct['Formula_ion'] = _usr_formula_charged
                            #             match_info_dct['Charge'] = _usr_charge
                            #             match_info_dct['MS1_obs_mz'] = _ms1_pr_mz
                            #             match_info_dct['MS1_obs_i'] = '%.2e' % float(_ms1_pr_i)
                            #             match_info_dct['Lib_mz'] = _usr_mz_lib
                            #             match_info_dct['MS2_scan_time'] = _usr_ms2_rt
                            #             match_info_dct['DDA#'] = _usr_ms2_dda_rank
                            #             match_info_dct['MS2_PR_mz'] = _usr_ms2_pr_mz
                            #             match_info_dct['Scan#'] = _usr_ms2_scan_id
                            #             match_info_dct['#Specific_peaks'] = (target_frag_count +
                            #                                                  target_nl_count)
                            #             match_info_dct['#Contaminated_peaks'] = (other_frag_count +
                            #                                                      other_nl_count)
                            #             match_info_dct['ppm'] = _exact_ppm
                            #             match_info_dct['SN_ratio'] = '%.1f' % usr_sn_ratio
                            #
                            #             # if any IO error while writing img output
                            #             if img_n == '-2':
                            #                 match_info_dct['img_name'] = '%s-2.%s' % (img_name_core[1:-4],
                            #                                                           img_typ)
                            #             else:
                            #                 match_info_dct['img_name'] = img_name_core[1:]
                            #
                            #             try:
                            #                 del match_info_dct['MATCH_INFO']
                            #                 del match_info_dct['OTHER_SIGNALS_INFO']
                            #                 del match_info_dct['MATCHED_FA_INFO']
                            #                 del match_info_dct['MATCHED_LYSO_INFO']
                            #             except KeyError:
                            #                 pass
                            #             _tmp_output_df = pd.DataFrame(data=match_info_dct, index=[tmp_idx])
                            #             tmp_df = tmp_df.append(_tmp_output_df)
                            #             tmp_idx += 1
    if tmp_df.shape[0] > 0:
        print('Size of the identified LPP_df %i, %i' % (tmp_df.shape[0], tmp_df.shape[1]))
        return tmp_df
    else:
        print('!! Size of the identified LPP_df == 0')
        return 'empty_df'
