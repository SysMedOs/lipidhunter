# -*- coding: utf-8 -*-
# Copyright 2015-2016 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de


import time
import pandas as pd

from LibLPPhunter.SpectraExtractor import extract_mzml
from LibLPPhunter.SpectraExtractor import get_spectra
from LibLPPhunter.SpectraExtractor import get_xic_all
from LibLPPhunter.ScoreGenerator import ScoreGenerator
from LibLPPhunter.PanelPloter import plot_spectra

start_time = time.clock()

# usr_lipid_type = 'PE'
# charge_mode = '[M-H]-'
usr_lipid_type = 'PC'
charge_mode = '[M+HCOO]-'
usr_mzml = r'D:\project_mzML\CM_DDA_neg_mzML\070120_CM_neg_70min_SIN_I.mzML'
usr_xlsx = r'D:\project_mzML\CM_DDA_neg_mzML\extractor_output\%s\%s_70min_SIN_I.xlsx' % (usr_lipid_type,
                                                                                         usr_lipid_type
                                                                                         )
output_folder = r'D:\project_mzML\CM_DDA_neg_mzML\images\%s\70min_SIN_I' % usr_lipid_type
output_sum_xlsx = r'D:\project_mzML\CM_DDA_neg_mzML\images\%s\70min_SIN_I\sum_%s_70min_SIN_I.xlsx' % (usr_lipid_type,
                                                                                                      usr_lipid_type
                                                                                                      )
fa_list_csv = r'D:\LPPhunter\FA_list.csv'
score_cfg = r'D:\LPPhunter\Score_cfg.xlsx'
key_frag_cfg = r'D:\LPPhunter\PL_specific_ion_cfg.xlsx'

usr_rt_range = [25, 25.3]
usr_pr_mz_range = [600, 1000]
usr_dda_top = 12
usr_ms1_precision = 50e-6
usr_ms2_precision = 500e-6

output_df = pd.DataFrame()

print('=== ==> --> Start to process')
print('=== ==> --> Phospholipid class: %s' % usr_lipid_type)

usr_df = pd.read_excel(usr_xlsx)
usr_df = usr_df.round({'mz': 6})
usr_df = usr_df.query('%f<= rt <= %f' % (usr_rt_range[0], usr_rt_range[1]))
usr_df = usr_df.query('%f<= mz <= %f' % (usr_pr_mz_range[0], usr_pr_mz_range[1]))
usr_df = usr_df.sort_values(by=['rt'])

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
                                                ms1_precision=usr_ms1_precision, msn_precision=usr_ms2_precision
                                                )

# remove bad precursors
checked_info_df = pd.DataFrame()
for _idx, _check_scan_se in usr_scan_info_df.iterrows():
    _function = _check_scan_se['function']
    _scan_id = _check_scan_se['scan_id']
    _tmp_usr_df = usr_df.query('function == %f and scan_id == %f' % (_function, _scan_id))
    checked_info_df = checked_info_df.append(_tmp_usr_df)

print('=== ==> --> Start to extract XIC')
xic_dct = get_xic_all(usr_df, usr_mzml, usr_rt_range, ms1_precision=usr_ms1_precision, msn_precision=500e-6)

print('=== ==> --> Number of XIC extracted: %i' % len(xic_dct.keys()))

# get spectra of one ABBR and plot
for _i, _row_se in checked_info_df.iterrows():
    print(_row_se)
    _usr_mz = _row_se['mz']
    _usr_ms1_obs_mz = _row_se['MS1_obs_mz']
    _usr_ms2_rt = _row_se['rt']
    _usr_abbr_bulk = _row_se['Abbreviation']
    _usr_formula = _row_se['Formula']
    _usr_ms2_function = _row_se['function']
    _usr_ms2_scan_id = _row_se['scan_id']
    _usr_mz_lib = _row_se['Lib_mz']
    _tmp_chk_df = usr_scan_info_df.query('pr_mz == %f and function == %i and scan_id == %i'
                                         % (_usr_mz, _usr_ms2_function, _usr_ms2_scan_id))

    _precision_checker = abs((_usr_ms1_obs_mz - _usr_mz_lib) / _usr_mz_lib)

    if _tmp_chk_df.shape[0] > 0 and _precision_checker <= usr_ms1_precision:
        _ms1_rt, _ms2_rt, ms1_spec_idx, ms2_spec_idx, ms1_df, ms2_df = get_spectra(_usr_mz,
                                                                                   _usr_ms2_function,
                                                                                   _usr_ms2_scan_id,
                                                                                   usr_scan_info_df,
                                                                                   usr_spectra_pl,
                                                                                   dda_top=usr_dda_top
                                                                                   )
        if ms1_df.shape[0] > 0 and ms2_df.shape[0] > 0:
            match_info_dct = score_calc.get_match(_usr_abbr_bulk, charge_mode, _usr_ms1_obs_mz, ms2_df,
                                                  ms2_precision=50e-6, ms2_threshold=100
                                                  )
            match_factor = match_info_dct['MATCH_INFO']
            score_df = match_info_dct['SCORE_INFO']
            fa_ident_df = match_info_dct['FA_INFO']
            lyso_ident_df = match_info_dct['LYSO_INFO']
            lyso_w_ident_df = match_info_dct['LYSO_W_INFO']

            if match_factor > 0 and score_df.shape[0] > 0:
                specific_check_dct = score_calc.get_specific_peaks(_usr_mz_lib, ms2_df,
                                                                   ms2_precision=200e-6, ms2_threshold=100
                                                                   )

                _usr_abbr_bulk = _usr_abbr_bulk.replace('(', '[')
                _usr_abbr_bulk = _usr_abbr_bulk.replace(')', ']')
                _usr_abbr_bulk = _usr_abbr_bulk.replace(':', '-')
                _usr_abbr_bulk = _usr_abbr_bulk.replace('\\', '_')
                _usr_abbr_bulk = _usr_abbr_bulk.replace('/', '_')

                print ('Check now for %s' % _usr_abbr_bulk)
                if fa_ident_df.shape[0]:
                    score_df = score_df[['Lipid_species', 'Score']]
                    score_df = score_df.rename({'Lipid_species': 'Proposed structures'})
                    score_df = score_df.query('Score >= 25')
                    score_df = score_df.sort_values(by='Score', ascending=False)
                    score_df = score_df.reset_index(drop=True)
                    score_df.index += 1
                    print(score_df)

                    # format fa info DataFrame
                    fa_ident_df = fa_ident_df[['Lipid_species', 'mz', 'i', 'ppm']].reset_index(drop=True)
                    fa_ident_df = fa_ident_df.rename({'Lipid_species': 'Identified species'})
                    fa_ident_df = fa_ident_df.round({'mz': 4, 'ppm': 2})
                    _fa_i_lst = []
                    for _idx, _fa_se in fa_ident_df.iterrows():
                        _fa_i_lst.append('%.2e' % float(_fa_se['i']))
                    fa_ident_df.loc[:, 'i'] = _fa_i_lst
                    fa_ident_df.index += 1
                    print(fa_ident_df)

                    # merge Lyso and Lyso - H2O
                    lyso_ident_df = lyso_ident_df.append(lyso_w_ident_df)
                    if lyso_ident_df.shape[0] > 0:
                        lyso_ident_df = lyso_ident_df.sort_values(by='i', ascending=False)
                        lyso_ident_df = lyso_ident_df[['Lipid_species', 'mz', 'i', 'ppm']].reset_index(drop=True)
                        lyso_ident_df = lyso_ident_df.rename({'Lipid_species': 'Identified species'})
                        lyso_ident_df = lyso_ident_df.round({'mz': 4, 'ppm': 2})
                        _lyso_i_lst = []
                        for _idx, _lyso_se in lyso_ident_df.iterrows():
                            _lyso_i_lst.append('%.2e' % float(_lyso_se['i']))
                        lyso_ident_df.loc[:, 'i'] = _lyso_i_lst
                        lyso_ident_df.index += 1
                        print(lyso_ident_df)
                    else:
                        lyso_ident_df = pd.DataFrame()

                    usr_ident_info_dct = {'SCORE_INFO': score_df, 'FA_INFO': fa_ident_df, 'LYSO_INFO': lyso_ident_df}

                    if score_df.shape[0] > 0:
                        img_name = output_folder + '\%.4f_rt%.4f_DDAtop%.0f_scan%.0f_%s.png' % (_usr_mz,
                                                                                                _usr_ms2_rt,
                                                                                                _usr_ms2_function - 1,
                                                                                                _usr_ms2_scan_id,
                                                                                                _usr_abbr_bulk
                                                                                                )

                        # print('_row_se', type(_row_se))
                        # print('xic_dct', type(xic_dct))
                        # print('xic_dct', _usr_ms1_obs_mz, type(xic_dct[_usr_ms1_obs_mz]))
                        # print('xic_dct', xic_dct[_usr_ms1_obs_mz].shape)
                        # print('xic_dct', xic_dct[_usr_ms1_obs_mz])
                        #
                        # print('_ms1_rt', type(_ms2_rt))
                        # print('ms1_df', type(ms1_df), ms1_df.shape)
                        # print('ms2_df', type(ms2_df), ms2_df.shape)
                        # print('img_name', type(img_name))
                        # print('usr_ms1_precision', type(usr_ms1_precision))
                        # print('usr_ident_info_dct', usr_ident_info_dct.keys())
                        # for _k in specific_check_dct.keys():
                        #     print(_k, type(specific_check_dct[_k]), specific_check_dct[_k].shape)

                        _ms1_pr_i, _ppm, isotope_checker = plot_spectra(_row_se, xic_dct, usr_ident_info_dct,
                                                                        _ms1_rt, _ms2_rt, ms1_df, ms2_df,
                                                                        specific_check_dct,
                                                                        save_img_as=img_name,
                                                                        ms1_precision=usr_ms1_precision
                                                                        )

                        if _ms1_pr_i > 0 and isotope_checker == 0:
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

                            _tmp_output_df['Bulk_identification'] = _usr_abbr_bulk
                            _tmp_output_df['MS1_obs_mz'] = _usr_ms1_obs_mz
                            _tmp_output_df['MS1_obs_i'] = '%.2e' % float(_ms1_pr_i)
                            _tmp_output_df['Lib_mz'] = _usr_mz_lib
                            _tmp_output_df['MS2_rt'] = _usr_ms2_rt
                            _tmp_output_df['MS2_function'] = _usr_ms2_function
                            _tmp_output_df['scan_id'] = _usr_ms2_scan_id
                            _tmp_output_df['specific peaks'] = target_frag_count + target_nl_count
                            _tmp_output_df['contaminated peaks'] = other_frag_count + other_nl_count
                            _tmp_output_df['ppm'] = _ppm

                            output_df = output_df.append(_tmp_output_df)

    else:
        print('!!!!!!!!!!!! PR NOT in list !!!!!!!!!!!!')

    print('---------------------NEXT---------------------------')
print('=== ==> --> Generate the output table')

# output_df['ppm'] = 1e6 * (output_df['MS1_obs_mz'] - output_df['Lib_mz']) / output_df['Lib_mz']
output_df = output_df.round({'MS1_obs_mz': 4, 'Lib_mz': 4, 'ppm': 2, 'MS2_rt': 3})
output_df['Proposed structures'] = output_df['Lipid_species']
output_df = output_df[['Bulk_identification', 'Proposed structures', 'Score', 'specific peaks', 'contaminated peaks',
                       'Lib_mz', 'MS1_obs_mz', 'MS1_obs_i', 'ppm', 'MS2_rt', 'MS2_function', 'scan_id']]
output_df = output_df.sort_values(by=['MS1_obs_mz', 'MS2_rt', 'Score'], ascending=[True, True, False])
output_df = output_df.reset_index(drop=True)
output_df.index += 1
output_df.to_excel(output_sum_xlsx)
print(output_sum_xlsx)
print('=== ==> --> saved >>> >>> >>>')

tot_run_time = time.clock() - start_time

print('>>> >>> >>> FINISHED in %f sec <<< <<< <<<' % tot_run_time)
