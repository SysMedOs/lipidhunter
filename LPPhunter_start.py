# -*- coding: utf-8 -*-
# Copyright 2015-2016 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

from __future__ import division

import time
import pandas as pd

from LibLipidHunter.SpectraExtractor import extract_mzml
from LibLipidHunter.SpectraExtractor import get_spectra
from LibLipidHunter.SpectraExtractor import get_xic_all
from LibLipidHunter.ScoreGenerator import ScoreGenerator
from LibLipidHunter.PanelPloter import plot_spectra
from LibLipidHunter.ScoreFilter import check_peaks
from LibLipidHunter.IsotopeHunter import IsotopeHunter

start_time = time.clock()

usr_lipid_type = 'PE'
charge_mode = '[M-H]-'
# usr_lipid_type = 'PC'
# charge_mode = '[M+HCOO]-'
usr_mzml = r'D:\project_mzML\CM_DDA_neg_mzML\070120_CM_neg_70min_SIN_II.mzML'
usr_xlsx = r'D:\project_mzML\CM_DDA_neg_mzML\extractor_output\%s\%s_70min_SIN_II_2.xlsx' % (usr_lipid_type,
                                                                                            usr_lipid_type
                                                                                            )

output_folder = r'D:\project_mzML\CM_DDA_neg_mzML\images\%s\70min_SIN_II_2' % usr_lipid_type
output_sum_xlsx = r'D:\project_mzML\CM_DDA_neg_mzML\images\%s\70min_SIN_II_2\sum_%s_70min_SIN_II.xlsx' % (
    usr_lipid_type,
    usr_lipid_type
)
fa_list_csv = r'D:\LPPhunter\FA_list.csv'
score_cfg = r'D:\LPPhunter\Score_cfg.xlsx'
key_frag_cfg = r'D:\LPPhunter\PL_specific_ion_cfg.xlsx'

usr_rt_range = [25, 26]
usr_pr_mz_range = [600, 1000]
usr_dda_top = 12
usr_ms1_threshold = 2000
usr_ms2_threshold = 100
usr_ms2_specific_peaks_threshold = 100
usr_ms1_precision = 50e-6
usr_ms2_precision = 200e-6
usr_ms2_specific_peaks_precision = 200e-6
usr_score_filter = 27.5
usr_isotope_score_filter = 85  # max 100

usr_ms1_ppm = int(usr_ms1_precision * 1e6)

output_df = pd.DataFrame()

print('=== ==> --> Start to process')
print('=== ==> --> Phospholipid class: %s' % usr_lipid_type)

usr_df = pd.read_excel(usr_xlsx)
usr_df = usr_df.round({'mz': 6})
usr_df = usr_df.round({'MS1_obs_mz': 6})
usr_df = usr_df.query('%f<= rt <= %f' % (usr_rt_range[0], usr_rt_range[1]))
usr_df = usr_df.query('%f<= mz <= %f' % (usr_pr_mz_range[0], usr_pr_mz_range[1]))
print(usr_df.shape)
usr_df['MS1_precision'] = (usr_df['MS1_obs_mz'] - usr_df['Lib_mz']) / usr_df['Lib_mz']
usr_df['ppm'] = 1e6 * usr_df['MS1_precision']
usr_df['abs_ppm'] = usr_df['ppm'].abs()
usr_df = usr_df.query('abs_ppm <= %i' % usr_ms1_ppm)
# usr_df = usr_df.sort_values(by=['Lib_mz', 'abs_ppm'], ascending=[True, True])
# usr_df = usr_df.drop_duplicates(subset=['Lib_mz', 'rt', 'function', 'scan_id'], keep='first')
usr_df = usr_df.sort_values(by=['rt'])
usr_df = usr_df.reset_index(drop=True)
print(usr_df.shape)

print('=== ==> --> Total precursor number: %i' % usr_df.shape[0])

# generate the indicator table

usr_fa_def_df = pd.read_csv(fa_list_csv)
usr_fa_def_df['C'] = usr_fa_def_df['C'].astype(int)
usr_fa_def_df['DB'] = usr_fa_def_df['DB'].astype(int)

usr_weight_df = pd.read_excel(score_cfg)

usr_key_frag_df = pd.read_excel(key_frag_cfg)
usr_key_frag_df = usr_key_frag_df.query('EXACTMASS > 0')
usr_key_frag_df = usr_key_frag_df[['CLASS', 'TYPE', 'EXACTMASS', 'PR_CHARGE', 'LABEL']]

score_calc = ScoreGenerator(usr_fa_def_df, usr_weight_df, usr_key_frag_df, usr_lipid_type)

print('=== ==> --> Start to parse mzML')
# extract all spectra from mzML to pandas DataFrame
usr_scan_info_df, usr_spectra_pl = extract_mzml(usr_mzml, usr_rt_range, dda_top=usr_dda_top,
                                                ms1_threshold=usr_ms1_threshold, ms2_threshold=usr_ms2_threshold,
                                                ms1_precision=usr_ms1_precision, ms2_precision=usr_ms2_precision
                                                )

usr_scan_info_df.to_excel('scan_info_df.xlsx')
# remove bad precursors
checked_info_df = pd.DataFrame()
for _idx, _check_scan_se in usr_scan_info_df.iterrows():
    _function = _check_scan_se['function']
    _scan_id = _check_scan_se['scan_id']
    _tmp_usr_df = usr_df.query('function == %f and scan_id == %f' % (_function, _scan_id))
    checked_info_df = checked_info_df.append(_tmp_usr_df)

checked_info_df.to_excel('checked_info_df.xlsx')
ms1_obs_mz_lst = usr_df['MS1_obs_mz'].tolist()
ms1_obs_mz_lst = set(ms1_obs_mz_lst)

print('=== ==> --> Start to extract XIC')
xic_dct = get_xic_all(usr_df, usr_mzml, usr_rt_range, ms1_precision=usr_ms1_precision, msn_precision=usr_ms2_precision)

print('=== ==> --> Number of XIC extracted: %i' % len(xic_dct.keys()))

plot_info_dct = {}
ms1_pr_mz_lst = []

# get spectra of one ABBR and plot
for _n, _subgroup_df in checked_info_df.groupby(['mz', 'Lib_mz', 'Formula', 'rt', 'Abbreviation']):
    _row_se = _subgroup_df.iloc[0, :]
    _usr_ms2_pr_mz = _row_se['mz']
    # _usr_ms1_obs_mz = _row_se['MS1_obs_mz']
    _usr_ms2_rt = _row_se['rt']
    _usr_formula = _row_se['Formula']
    _usr_ms2_function = _row_se['function']
    _usr_ms2_scan_id = _row_se['scan_id']
    _usr_mz_lib = _row_se['Lib_mz']
    _usr_abbr_bulk = _row_se['Abbreviation']
    _tmp_chk_df = usr_scan_info_df.query('pr_mz == %.6f and function == %i and scan_id == %i'
                                         % (_usr_ms2_pr_mz, _usr_ms2_function, _usr_ms2_scan_id))

    if _tmp_chk_df.shape[0] == 1:
        print('>>> >>> >>> Processing:', _tmp_chk_df.head())
        print('>>> >>> >>> >>> MS2 PR m/z %f' % _usr_ms2_pr_mz)
        spec_info_dct = get_spectra(_usr_ms2_pr_mz, _usr_mz_lib, _usr_ms2_function, _usr_ms2_scan_id, ms1_obs_mz_lst,
                                    usr_scan_info_df, usr_spectra_pl,
                                    dda_top=usr_dda_top, ms1_precision=usr_ms1_precision
                                    )
        _ms1_pr_i = spec_info_dct['ms1_i']
        _ms1_pr_mz = spec_info_dct['ms1_mz']
        _ms1_rt = spec_info_dct['ms1_rt']
        _ms2_rt = spec_info_dct['ms2_rt']
        _ms1_spec_idx = spec_info_dct['_ms1_spec_idx']
        _ms2_spec_idx = spec_info_dct['_ms2_spec_idx']
        _ms1_df = spec_info_dct['_ms1_df']
        _ms2_df = spec_info_dct['_ms2_df']

        # _ms1_subgroup_df = _subgroup_df.query('MS1_obs_mz == %f' % _ms1_pr_mz)
        # if _ms1_subgroup_df.shape[0] > 0 and _ms1_df.shape[0] > 0 and _ms2_df.shape[0] > 0:
        if _ms1_df.shape[0] > 0 and _ms2_df.shape[0] > 0:

            print('>>> >>> >>> >>> Best PR on MS1: %f' % _ms1_pr_mz)
            _row_se['MS1_obs_mz'] = _ms1_pr_mz
            print('>>> >>> >>> >>> Entry Info >>> >>> >>> >>> ')
            print(_row_se)
            match_info_dct = score_calc.get_match(_usr_abbr_bulk, charge_mode, _ms1_pr_mz, _ms2_df,
                                                  ms2_precision=usr_ms2_precision, ms2_threshold=usr_ms2_threshold
                                                  )
            match_factor = match_info_dct['MATCH_INFO']
            score_df = match_info_dct['SCORE_INFO']
            fa_ident_df = match_info_dct['FA_INFO']
            lyso_ident_df = match_info_dct['LYSO_INFO']
            lyso_w_ident_df = match_info_dct['LYSO_W_INFO']

            if match_factor > 0 and score_df.shape[0] > 0 and fa_ident_df.shape[0]:
                specific_check_dct = score_calc.get_specific_peaks(_usr_mz_lib, _ms2_df,
                                                                   ms2_precision=usr_ms2_specific_peaks_precision,
                                                                   ms2_threshold=usr_ms2_specific_peaks_threshold
                                                                   )

                _usr_abbr_bulk = _usr_abbr_bulk.replace('(', '[')
                _usr_abbr_bulk = _usr_abbr_bulk.replace(')', ']')
                _usr_abbr_bulk = _usr_abbr_bulk.replace(':', '-')
                _usr_abbr_bulk = _usr_abbr_bulk.replace('\\', '_')
                _usr_abbr_bulk = _usr_abbr_bulk.replace('/', '_')

                print ('>>> >>> Check now for bulk identification as %s' % _usr_abbr_bulk)

                usr_ident_info_dct = check_peaks(score_df, fa_ident_df, lyso_ident_df, lyso_w_ident_df,
                                                 score_filter=usr_score_filter)

                score_df = usr_ident_info_dct['SCORE_INFO']
                if score_df.shape[0] > 0 and _ms1_pr_i > 0:
                    img_name = output_folder + '\%.4f_rt%.4f_DDAtop%.0f_scan%.0f_%s.png' % (_usr_ms2_pr_mz,
                                                                                            _usr_ms2_rt,
                                                                                            _usr_ms2_function - 1,
                                                                                            _usr_ms2_scan_id,
                                                                                            _usr_abbr_bulk
                                                                                            )

                    isotope_hunter = IsotopeHunter()
                    isotope_score, isotope_checker_dct = isotope_hunter.get_isotope_score(_ms1_pr_mz, _ms1_pr_i,
                                                                                          _usr_formula, _ms1_df,
                                                                                          isotope_number=2)
                    print('isotope_score: %f' % isotope_score)
                    if isotope_score >= usr_isotope_score_filter:
                        _ms1_pr_i, _ppm, isotope_checker, isotope_score = plot_spectra(_row_se, xic_dct,
                                                                                       usr_ident_info_dct,
                                                                                       _ms1_rt, _ms2_rt, _ms1_df,
                                                                                       _ms2_df,
                                                                                       specific_check_dct,
                                                                                       isotope_checker_dct,
                                                                                       isotope_score,
                                                                                       save_img_as=img_name,
                                                                                       ms1_precision=usr_ms1_precision
                                                                                       )

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
                            else:
                                target_frag_count = 0
                            if 'TARGET_NL' in specific_check_dct.keys():
                                target_nl_df = specific_check_dct['TARGET_NL']
                                target_nl_count = target_nl_df.shape[0]
                            else:
                                target_nl_count = 0

                            _tmp_output_df['Bulk identification'] = _usr_abbr_bulk
                            _tmp_output_df['MS1_obs_mz'] = _ms1_pr_mz
                            _tmp_output_df['MS1_obs_i'] = '%.2e' % float(_ms1_pr_i)
                            _tmp_output_df['Lib_mz'] = _usr_mz_lib
                            _tmp_output_df['MS2_rt'] = _usr_ms2_rt
                            _tmp_output_df['MS2_function'] = _usr_ms2_function
                            _tmp_output_df['MS2_PR_MZ'] = _usr_ms2_pr_mz
                            _tmp_output_df['Scan#'] = _usr_ms2_scan_id
                            _tmp_output_df['Specific peaks'] = target_frag_count + target_nl_count
                            _tmp_output_df['Contaminated peaks'] = other_frag_count + other_nl_count
                            _tmp_output_df['ppm'] = _ppm
                            _tmp_output_df['Isotope_score'] = '%.2f' % isotope_score

                            output_df = output_df.append(_tmp_output_df)

print('=== ==> --> Generate the output table')

# output_df['ppm'] = 1e6 * (output_df['MS1_obs_mz'] - output_df['Lib_mz']) / output_df['Lib_mz']
output_df = output_df.round({'MS1_obs_mz': 4, 'Lib_mz': 4, 'ppm': 2, 'MS2_rt': 3})
output_df['Proposed structures'] = output_df['Lipid_species']
output_df = output_df[['Bulk identification', 'Proposed structures', 'Score', 'Specific peaks', 'Contaminated peaks',
                       'Lib_mz', 'MS1_obs_mz', 'MS1_obs_i', 'ppm', 'Isotope_score',
                       'MS2_PR_MZ', 'MS2_rt', 'MS2_function', 'Scan#']]
output_df = output_df.sort_values(by=['MS1_obs_mz', 'MS2_rt', 'Score'], ascending=[True, True, False])
output_df = output_df.reset_index(drop=True)
output_df.index += 1
output_df.to_excel(output_sum_xlsx)
print(output_sum_xlsx)
print('=== ==> --> saved >>> >>> >>>')

tot_run_time = time.clock() - start_time

print('>>> >>> >>> FINISHED in %f sec <<< <<< <<<' % tot_run_time)
