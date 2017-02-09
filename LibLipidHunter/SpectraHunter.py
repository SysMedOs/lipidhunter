# -*- coding: utf-8 -*-
# Copyright 2016-2017 LPP team, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of LipidHunter.
# For more info please contact:
#     LPP team oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
#     Developer Georgia Angelidou georgia.angelidou@uni-leipzig.de

from __future__ import division

import time
import os
import pandas as pd

from LibLipidHunter.SpectraExtractor import extract_mzml
from LibLipidHunter.SpectraExtractor import get_spectra
from LibLipidHunter.SpectraExtractor import get_xic_all
from LibLipidHunter.ScoreGenerator import ScoreGenerator
from LibLipidHunter.PanelPloter import plot_spectra
from LibLipidHunter.ScoreFilter import check_peaks
from LibLipidHunter.IsotopeHunter import IsotopeHunter
from LibLipidHunter.AbbrElemCalc import BulkAbbrFormula
from LibLipidHunter.LogPageCreator import LogPageCreator


def huntlipids(param_dct):
    """

    :param param_dct:
    :return:
    """

    start_time = time.clock()

    usr_lipid_type = param_dct['lipid_type']
    charge_mode = param_dct['charge_mode']
    usr_vendor = param_dct['vendor']
    usr_xlsx = param_dct['lipids_info_path_str']
    usr_mzml = param_dct['mzml_path_str']
    output_folder = param_dct['img_output_folder_str']
    output_sum_xlsx = param_dct['xlsx_output_path_str']
    fa_list_csv = param_dct['fa_white_list_cfg']

    key_frag_cfg = param_dct['lipid_specific_cfg']
    score_cfg = param_dct['score_cfg']

    usr_rt_range = [param_dct['rt_start'], param_dct['rt_end']]
    usr_pr_mz_range = [param_dct['mz_start'], param_dct['mz_end']]
    usr_dda_top = param_dct['dda_top']
    usr_ms1_threshold = param_dct['ms_th']
    usr_ms2_threshold = param_dct['ms2_th']
    usr_ms2_hg_threshold = param_dct['hg_th']
    usr_ms1_precision = param_dct['ms_ppm'] * 1e-6
    usr_ms2_precision = param_dct['ms2_ppm'] * 1e-6
    usr_ms2_hg_precision = param_dct['hg_ppm'] * 1e-6
    usr_score_filter = param_dct['score_filter']
    usr_isotope_score_filter = param_dct['isotope_score_filter']
    usr_ms2_info_th = param_dct['ms2_infopeak_threshold']
    usr_ms2_hginfo_th = param_dct['ms2_hginfopeak_threshold']
    usr_rank_mode = param_dct['rank_score']
    usr_fast_isotope = param_dct['fast_isotope']

    # hunter_folder = param_dct['hunter_folder']

    if usr_rank_mode is True:
        score_mode = 'Rank mode'
    else:
        score_mode = 'Relative intensity mode'

    if usr_fast_isotope is True:
        isotope_score_mode = '(Fast mode)'
    else:
        isotope_score_mode = ''

    usr_ms1_ppm = int(param_dct['ms_ppm'])

    hunter_start_time_str = param_dct['hunter_start_time']
    isotope_hunter = IsotopeHunter()
    abbr2formula = BulkAbbrFormula()

    # keep stay in current working directory
    current_path = os.getcwd()
    if os.path.isdir(output_folder):
        os.chdir(output_folder)
        if os.path.isdir('LipidHunter_Results_Figures_%s' % hunter_start_time_str):
            pass
        else:
            os.mkdir('LipidHunter_Results_Figures_%s' % hunter_start_time_str)
    os.chdir(current_path)

    log_pager = LogPageCreator(output_folder, hunter_start_time_str, param_dct)

    output_df = pd.DataFrame()

    print('=== ==> --> Start to process')
    print('=== ==> --> Phospholipid class: %s' % usr_lipid_type)

    usr_df = pd.read_excel(usr_xlsx)
    # usr_df = usr_df.round({'mz': 6})
    # usr_df = usr_df.round({'MS1_obs_mz': 6})
    tmp_usr_df = usr_df.copy()
    tmp_usr_df = tmp_usr_df.query('%f<= scan_time <= %f' % (usr_rt_range[0], usr_rt_range[1]))
    tmp_usr_df = tmp_usr_df.query('%f<= MS2_PR_mz <= %f' % (usr_pr_mz_range[0], usr_pr_mz_range[1]))
    usr_df = tmp_usr_df.copy()

    for _i, _r in usr_df.iterrows():
        usr_df.set_value(_i, 'MS2_PR_mz', round(_r['MS2_PR_mz'], 6))
        usr_df.set_value(_i, 'MS1_obs_mz', round(_r['MS1_obs_mz'], 6))
        usr_df.set_value(_i, 'MS1_precision', (_r['MS1_obs_mz'] - _r['Lib_mz']) / _r['Lib_mz'])
        usr_df.set_value(_i, 'ppm', 1e6 * (_r['MS1_obs_mz'] - _r['Lib_mz']) / _r['Lib_mz'])
        usr_df.set_value(_i, 'abs_ppm', abs(1e6 * (_r['MS1_obs_mz'] - _r['Lib_mz']) / _r['Lib_mz']))
    tmp_usr_df = usr_df.query('abs_ppm <= %i' % usr_ms1_ppm)
    # usr_df = usr_df.sort_values(by=['Lib_mz', 'abs_ppm'], ascending=[True, True])
    # usr_df = usr_df.drop_duplicates(subset=['Lib_mz', 'rt', 'DDA_rank', 'scan_id'], keep='first')
    usr_df = tmp_usr_df.sort_values(by=['scan_time'])
    tmp_usr_df = usr_df.reset_index(drop=True)
    usr_df = tmp_usr_df.copy()

    print('=== ==> --> Total precursor number: %i' % usr_df.shape[0])

    # generate the indicator table

    usr_fa_def_df = pd.read_csv(fa_list_csv)
    usr_fa_def_df.loc[:, 'C'] = usr_fa_def_df['C'].astype(int)
    usr_fa_def_df.loc[:, 'DB'] = usr_fa_def_df['DB'].astype(int)

    usr_weight_df = pd.read_excel(score_cfg)

    usr_key_frag_df = pd.read_excel(key_frag_cfg)
    usr_key_frag_df = usr_key_frag_df.query('EXACTMASS > 0')

    # get the information from the following columns and leave the rewark back
    usr_key_frag_df = usr_key_frag_df[['CLASS', 'TYPE', 'EXACTMASS', 'PR_CHARGE', 'LABEL', 'CHARGE_MODE']]

    score_calc = ScoreGenerator(usr_fa_def_df, usr_weight_df, usr_key_frag_df, usr_lipid_type, ion_charge=charge_mode)

    print('=== ==> --> Start to parse mzML')
    # extract all spectra from mzML to pandas DataFrame
    usr_scan_info_df, usr_spectra_pl = extract_mzml(usr_mzml, usr_rt_range, dda_top=usr_dda_top,
                                                    ms1_threshold=usr_ms1_threshold, ms2_threshold=usr_ms2_threshold,
                                                    ms1_precision=usr_ms1_precision, ms2_precision=usr_ms2_precision,
                                                    vendor=usr_vendor
                                                    )

    # remove bad precursors
    checked_info_df = pd.DataFrame()
    for _idx, _check_scan_se in usr_scan_info_df.iterrows():
        _dda_rank = _check_scan_se['DDA_rank']
        _scan_id = _check_scan_se['scan_number']
        _tmp_usr_df = usr_df.query('DDA_rank == %f and scan_number == %f' % (_dda_rank, _scan_id))
        checked_info_df = checked_info_df.append(_tmp_usr_df)

    ms1_obs_mz_lst = usr_df['MS1_obs_mz'].tolist()
    ms1_obs_mz_lst = set(ms1_obs_mz_lst)

    print('=== ==> --> Start to extract XIC')
    try:
        xic_dct = get_xic_all(usr_df, usr_mzml, usr_rt_range, ms1_precision=usr_ms1_precision,
                              msn_precision=usr_ms2_precision, vendor=usr_vendor)

    except KeyError:
        return u'Nothing found! Check mzML vendor settings!'

    print('=== ==> --> Number of XIC extracted: %i' % len(xic_dct.keys()))

    # plot_info_dct = {}
    # ms1_pr_mz_lst = []
    target_ident_lst = []
    ident_page_idx = 1

    # get spectra of one ABBR and plot
    for _n, _subgroup_df in checked_info_df.groupby(['MS2_PR_mz', 'Lib_mz', 'Formula', 'scan_time', 'Abbreviation']):
        _row_se = _subgroup_df.iloc[0, :].squeeze()
        _usr_ms2_pr_mz = _row_se['MS2_PR_mz']
        # _usr_ms1_obs_mz = _row_se['MS1_obs_mz']
        _usr_ms2_rt = _row_se['scan_time']
        # _usr_formula_charged = _row_se['Formula']
        _usr_charge = _row_se['Ion']
        _usr_ms2_dda_rank = _row_se['DDA_rank']
        _usr_ms2_scan_id = _row_se['scan_number']
        _usr_mz_lib = _row_se['Lib_mz']
        _usr_abbr_bulk = _row_se['Abbreviation']
        _tmp_chk_df = usr_scan_info_df.query('MS2_PR_mz == %.6f and DDA_rank == %i and scan_number == %i'
                                             % (_usr_ms2_pr_mz, _usr_ms2_dda_rank, _usr_ms2_scan_id))

        _usr_formula_charged, usr_elem_charged_dct = abbr2formula.get_formula(_usr_abbr_bulk, charge=_usr_charge)
        _usr_formula, usr_elem_dct = abbr2formula.get_formula(_usr_abbr_bulk, charge='')

        if _tmp_chk_df.shape[0] == 1:
            print('>>> >>> >>> Processing:', _tmp_chk_df.head())
            print('>>> >>> >>> >>> MS2 PR m/z %f' % _usr_ms2_pr_mz)
            usr_spec_info_dct = get_spectra(_usr_ms2_pr_mz, _usr_mz_lib, _usr_ms2_dda_rank, _usr_ms2_scan_id,
                                            ms1_obs_mz_lst, usr_scan_info_df, usr_spectra_pl,
                                            dda_top=usr_dda_top, ms1_precision=usr_ms1_precision
                                            )
            _ms1_pr_i = usr_spec_info_dct['ms1_i']
            _ms1_pr_mz = usr_spec_info_dct['ms1_mz']
            _ms1_df = usr_spec_info_dct['ms1_df']
            _ms2_df = usr_spec_info_dct['ms2_df']
            # _ms1_rt = usr_spec_info_dct['ms1_rt']
            # _ms2_rt = usr_spec_info_dct['ms2_rt']
            # _ms1_pr_df = usr_spec_info_dct['ms1_pr_df']
            # _ms1_spec_idx = usr_spec_info_dct['_ms1_spec_idx']
            # _ms2_spec_idx = usr_spec_info_dct['_ms2_spec_idx']

            if _ms1_pr_mz > 0.0 and _ms1_df.shape[0] > 0 and _ms2_df.shape[0] > 0:

                print('>>> >>> >>> >>> Best PR on MS1: %f' % _ms1_pr_mz)

                isotope_score, isotope_checker_dct = isotope_hunter.get_isotope_score(_ms1_pr_mz, _ms1_pr_i,
                                                                                      _usr_formula_charged, _ms1_df,
                                                                                      isotope_number=2,
                                                                                      only_c=usr_fast_isotope)
                print('isotope_score: %f' % isotope_score)
                if isotope_score >= usr_isotope_score_filter:
                    print('>>> isotope_check PASSED! >>> >>> >>>')
                    print('>>> >>> >>> >>> Entry Info >>> >>> >>> >>> ')
                    _row_se.set_value('MS1_obs_mz', _ms1_pr_mz)
                    _exact_ppm = 1e6 * (_ms1_pr_mz - _usr_mz_lib) / _usr_mz_lib
                    _row_se.set_value('ppm', _exact_ppm)
                    _row_se.set_value('abs_ppm', abs(_exact_ppm))
                    print(_row_se)
                    match_info_dct = score_calc.get_match(_usr_abbr_bulk, charge_mode, _ms1_pr_mz, _ms2_df,
                                                          ms2_precision=usr_ms2_precision,
                                                          ms2_threshold=usr_ms2_threshold,
                                                          ms2_infopeak_threshold=usr_ms2_info_th,
                                                          rank_mode=usr_rank_mode
                                                          )
                    match_factor = match_info_dct['MATCH_INFO']
                    score_df = match_info_dct['SCORE_INFO']
                    fa_ident_df = match_info_dct['FA_INFO']
                    lyso_ident_df = match_info_dct['LYSO_INFO']
                    lyso_w_ident_df = match_info_dct['LYSO_W_INFO']

                    if match_factor > 0 and score_df.shape[0] > 0 and fa_ident_df.shape[0] > 0:

                        usr_ident_info_dct = check_peaks(score_df, fa_ident_df, lyso_ident_df, lyso_w_ident_df,
                                                         score_filter=usr_score_filter)

                        score_df = usr_ident_info_dct['SCORE_INFO']
                        if score_df.shape[0] > 0 and _ms1_pr_i > 0:
                            print ('>>> >>> Check now for bulk identification as %s' % _usr_abbr_bulk)

                            specific_check_dct = score_calc.get_specific_peaks(_usr_mz_lib, _ms2_df,
                                                                               ms2_precision=usr_ms2_hg_precision,
                                                                               ms2_threshold=usr_ms2_hg_threshold,
                                                                               ms2_hginfo_threshold=usr_ms2_hginfo_th,
                                                                               vendor=usr_vendor
                                                                               )
                            # format abbr. for file names
                            _save_abbr_bulk = _usr_abbr_bulk
                            _save_abbr_bulk = _save_abbr_bulk.replace('(', '[')
                            _save_abbr_bulk = _save_abbr_bulk.replace(')', ']')
                            _save_abbr_bulk = _save_abbr_bulk.replace(':', '-')
                            _save_abbr_bulk = _save_abbr_bulk.replace('\\', '_')
                            _save_abbr_bulk = _save_abbr_bulk.replace('/', '_')

                            img_name_core = ('\%.4f_rt%.3f_DDAtop%.0f_scan%.0f_%s.png'
                                             % (_usr_ms2_pr_mz, _usr_ms2_rt, _usr_ms2_dda_rank,
                                                _usr_ms2_scan_id, _save_abbr_bulk)
                                             )

                            img_name = (output_folder + r'\LipidHunter_Results_Figures_%s' % hunter_start_time_str
                                        + img_name_core)

                            isotope_checker, isotope_score = plot_spectra(_row_se, xic_dct, usr_ident_info_dct,
                                                                          usr_spec_info_dct, specific_check_dct,
                                                                          isotope_checker_dct, isotope_score,
                                                                          _usr_formula_charged, _usr_charge,
                                                                          save_img_as=img_name,
                                                                          ms1_precision=usr_ms1_precision,
                                                                          score_mode=score_mode,
                                                                          isotope_mode=isotope_score_mode
                                                                          )

                            print('==> check for output -->')

                            if _ms1_pr_i > 0 and isotope_checker == 0 and isotope_score > usr_isotope_score_filter:
                                _tmp_output_df = score_df

                                if 'OTHER_FRAG' in specific_check_dct.keys():
                                    other_frag_df = specific_check_dct['OTHER_FRAG']
                                    other_frag_count = other_frag_df.shape[0]
                                else:
                                    other_frag_count = 0
                                if 'OTHER_NL' in specific_check_dct.keys():
                                    other_nl_df = specific_check_dct['OTHER_NL']
                                    other_nl_count = other_nl_df.shape[0]
                                else:
                                    other_nl_count = 0
                                if 'TARGET_FRAG' in specific_check_dct.keys():
                                    target_frag_df = specific_check_dct['TARGET_FRAG']
                                    target_frag_count = target_frag_df.shape[0]
                                    target_frag_col_lst = target_frag_df.columns.tolist()
                                    for _frag_abbr in target_frag_col_lst:
                                        if _frag_abbr not in ['mz', 'i', 'LABEL', 'CLASS']:
                                            for _i, _f_se in target_frag_df.iterrows():
                                                if _f_se['LABEL'] == _frag_abbr:
                                                    _tmp_output_df[_frag_abbr] = _f_se[_frag_abbr]
                                                    if _frag_abbr not in target_ident_lst:
                                                        target_ident_lst.append(_frag_abbr)
                                else:
                                    target_frag_count = 0
                                if 'TARGET_NL' in specific_check_dct.keys():
                                    target_nl_df = specific_check_dct['TARGET_NL']
                                    target_nl_count = target_nl_df.shape[0]
                                    target_nl_col_lst = target_nl_df.columns.tolist()
                                    for _nl_abbr in target_nl_col_lst:
                                        if _nl_abbr not in ['mz', 'i', 'LABEL', 'CLASS']:
                                            for _i, _n_se in target_nl_df.iterrows():
                                                if _n_se['LABEL'] == _nl_abbr:
                                                    _tmp_output_df[_nl_abbr] = _n_se[_nl_abbr]
                                                    if _nl_abbr not in target_ident_lst:
                                                        target_ident_lst.append(_nl_abbr)
                                else:
                                    target_nl_count = 0

                                _tmp_output_df['Bulk_identification'] = _usr_abbr_bulk
                                _tmp_output_df['Formula_neutral'] = _usr_formula
                                _tmp_output_df['Formula_ion'] = _usr_formula_charged
                                _tmp_output_df['Charge'] = _usr_charge
                                _tmp_output_df['MS1_obs_mz'] = _ms1_pr_mz
                                _tmp_output_df['MS1_obs_i'] = '%.2e' % float(_ms1_pr_i)
                                _tmp_output_df['Lib_mz'] = _usr_mz_lib
                                _tmp_output_df['MS2_scan_time'] = _usr_ms2_rt
                                _tmp_output_df['DDA#'] = _usr_ms2_dda_rank - 1
                                _tmp_output_df['MS2_PR_mz'] = _usr_ms2_pr_mz
                                _tmp_output_df['Scan#'] = _usr_ms2_scan_id
                                _tmp_output_df['#Specific_peaks'] = target_frag_count + target_nl_count
                                _tmp_output_df['#Contaminated_peaks'] = other_frag_count + other_nl_count
                                _tmp_output_df['ppm'] = _exact_ppm
                                _tmp_output_df['Isotope_score'] = '%.2f' % isotope_score

                                output_df = output_df.append(_tmp_output_df)

                                log_pager.add_info(img_name_core, ident_page_idx, _tmp_output_df)
                                ident_page_idx += 1

    print('=== ==> --> Generate the output table')
    if output_df.shape[0] > 0:
        output_header_lst = output_df.columns.tolist()
        if usr_lipid_type == 'TG':
            for _i_check in ['i_sn1', 'i_sn2', 'i_sn3', 'i_M-sn1', 'i_M-sn2', 'i_M-sn3', 'i_M-(sn1+sn2)',
                             'i_M-(sn1+sn3)', 'i_M-(sn2+sn3)']:
                if _i_check not in output_header_lst:
                    output_df[_i_check] = 0.0
            output_round_dct = {r'MS1_obs_mz': 4, r'Lib_mz': 4, 'ppm': 2, 'MS2_scan_time': 3,
                                'i_sn1': 2, 'i_sn2': 2, 'i_sn3': 2, 'i_M-sn1': 2, 'i_M-sn2': 2,
                                'i_M-sn3': 2, 'i_M-(sn1+sn2)': 2, 'i_M-(sn2+sn3)': 2, 'i_M-(sn1+sn3)': 2
                                }
        else:
            for _i_check in ['i_sn1', 'i_sn2', 'i_M-sn1', 'i_M-sn2', 'i_M-(sn1-H2O)', 'i_M-(sn2-H2O)']:
                if _i_check not in output_header_lst:
                    output_df[_i_check] = 0.0

            output_round_dct = {r'MS1_obs_mz': 4, r'Lib_mz': 4, 'ppm': 2, 'MS2_scan_time': 3,
                                'i_sn1': 2, 'i_sn2': 2, 'i_M-sn1': 2, 'i_M-sn2': 2,
                                'i_M-(sn1-H2O)': 2, 'i_M-(sn2-H2O)': 2
                                }
        # add intensities of target peaks to round list
        if len(target_ident_lst) > 0:
            for _t in target_ident_lst:
                output_round_dct[_t] = 2
        output_df = output_df.round(output_round_dct)

        # output_df['Proposed structures'] = output_df['Lipid_species']
        print(output_df.head())
        output_df.rename(columns={'Score': 'LipidHunter_Score',
                                  '#Contaminated_peaks': '#Unspecific_peaks'}, inplace=True)
        if usr_lipid_type == 'TG':
            output_header_lst = ['Bulk_identification', 'Proposed_structures', 'Formula_neutral', 'Formula_ion',
                                 'Charge', 'Lib_mz', 'ppm', 'LipidHunter_Score', 'MS1_obs_mz', 'MS1_obs_i',
                                 'Isotope_score', r'MS2_PR_mz', 'MS2_scan_time', 'DDA#', 'Scan#',
                                 'i_sn1', 'i_sn2', 'i_sn3', 'i_M-sn1', 'i_M-sn2', 'i_M-sn3', 'i_M-(sn1+sn2)',
                                 'i_M-(sn2+sn3)', 'i_M-(sn1+sn3)', '#Specific_peaks']
        else:
            output_header_lst = ['Bulk_identification', 'Proposed_structures', 'Formula_neutral', 'Formula_ion',
                                 'Charge', 'Lib_mz', 'ppm', 'LipidHunter_Score', 'MS1_obs_mz', 'MS1_obs_i',
                                 'Isotope_score', r'MS2_PR_mz', 'MS2_scan_time', 'DDA#', 'Scan#', 'i_sn1', 'i_sn2',
                                 'i_M-sn1', 'i_M-sn2', 'i_M-(sn1-H2O)', 'i_M-(sn2-H2O)', '#Specific_peaks']
        output_header_lst += target_ident_lst
        output_header_lst += ['#Unspecific_peaks']
        print(output_df.head())
        print(output_header_lst)
        output_df = output_df[output_header_lst]
        output_df = output_df.sort_values(by=['MS1_obs_mz', 'MS2_scan_time', 'LipidHunter_Score'],
                                          ascending=[True, True, False])
        output_df = output_df.reset_index(drop=True)
        output_df.index += 1
        output_df.to_excel(output_sum_xlsx, index=False)
        print(output_sum_xlsx)
        print('=== ==> --> saved >>> >>> >>>')

    log_pager.close_page()

    tot_run_time = time.clock() - start_time

    print('>>> >>> >>> FINISHED in %f sec <<< <<< <<<' % tot_run_time)

    return tot_run_time


if __name__ == '__main__':
    pass
