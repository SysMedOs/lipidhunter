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

from LibLPPhunter.FAindicator import FAindicator
from LibLPPhunter.FAindicator import Lyso_indicator

st_time = time.clock()

pl_type = 'PC'
charge_mode = '[M+FA-H]-'
usr_mzml = r'D:\project_mzML\CM_DDA_neg_mzML\070120_CM_neg_70min_SIN_II.mzML'
usr_xlsx = r'D:\project_mzML\CM_DDA_neg_mzML\extractor_output\%s\%s_70min_SIN_II.xlsx' % (pl_type, pl_type)
output_folder = r'D:\project_mzML\CM_DDA_neg_mzML\images\%s\70min_SIN_II' % pl_type
fa_list_csv = r'D:\LPPhunter\FA_list.csv'
score_cfg = r'D:\LPPhunter\Score_cfg.xlsx'

usr_rt_range = [25, 26]
usr_pr_mz_range = [600, 1000]
usr_dda_top = 12
usr_ms1_precision = 50e-6
usr_ms2_precision = 500e-6

print('=== ==> --> Start to process')
print('=== ==> --> Phospholipid class: %s' % pl_type)

usr_df = pd.read_excel(usr_xlsx)
usr_df = usr_df.round({'mz': 6})
usr_df = usr_df.query('%f<= rt <= %f' % (usr_rt_range[0], usr_rt_range[1]))
usr_df = usr_df.query('%f<= mz <= %f' % (usr_pr_mz_range[0], usr_pr_mz_range[1]))
usr_df = usr_df.sort_values(by=['rt'])

print('=== ==> --> Total precursor number: %i' % usr_df.shape[0])

# generate the indicator table
fa_indicator = FAindicator(fa_list_csv)
lyso_indicator = Lyso_indicator(fa_list_csv)
usr_weight_df = pd.read_excel(score_cfg)
usr_fa_def_df = pd.read_csv(fa_list_csv)
usr_fa_def_df['C'] = usr_fa_def_df['C'].astype(int)
usr_fa_def_df['DB'] = usr_fa_def_df['DB'].astype(int)

score_calc = ScoreGenerator(usr_fa_def_df, usr_weight_df)

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
xic_dct = get_xic_all(usr_df, usr_mzml, usr_rt_range, ms1_precision=50e-6, msn_precision=500e-6)

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

    if _tmp_chk_df.shape[0] > 0:
        _ms1_rt, _ms2_rt, ms1_spec_idx, ms2_spec_idx, ms1_df, ms2_df = get_spectra(_usr_mz,
                                                                                   _usr_ms2_function,
                                                                                   _usr_ms2_scan_id,
                                                                                   usr_scan_info_df,
                                                                                   usr_spectra_pl,
                                                                                   dda_top=usr_dda_top
                                                                                   )

        score_df, fa_ident_df, lyso_ident_df, lyso_w_ident_df = score_calc.get_match(_usr_abbr_bulk, charge_mode,
                                                                                     _usr_ms1_obs_mz, ms2_df,
                                                                                     ms2_precision=50e-6,
                                                                                     ms2_threshold=100
                                                                                     )

        _usr_abbr_bulk = _usr_abbr_bulk.replace('(', '[')
        _usr_abbr_bulk = _usr_abbr_bulk.replace(')', ']')
        _usr_abbr_bulk = _usr_abbr_bulk.replace(':', '-')
        _usr_abbr_bulk = _usr_abbr_bulk.replace('\\', '_')
        _usr_abbr_bulk = _usr_abbr_bulk.replace('/', '_')

        print ('_usr_abbr_bulk', _usr_abbr_bulk)
        score_df = score_df[['Lipid_abbr', 'Score']]
        score_df = score_df.rename({'Lipid_abbr': 'Proposed structure'})
        score_df = score_df.query('Score > 20')
        score_df = score_df.sort_values(by='Score', ascending=False)
        score_df = score_df.reset_index(drop=True)
        print('score_df')
        print(score_df)
        lyso_ident_df = lyso_ident_df.append(lyso_w_ident_df)
        lyso_ident_df = lyso_ident_df.sort_values(by='i', ascending=False)

        lyso_ident_df = lyso_ident_df.reset_index(drop=True)
        fa_ident_df = fa_ident_df.reset_index(drop=True)

        lyso_ident_df = lyso_ident_df.round({'mz': 4, 'ppm': 2})
        fa_ident_df = fa_ident_df.round({'mz': 4, 'ppm': 2})
        usr_ident_info_df = {'SCORE_INFO': score_df, 'FA_INFO': fa_ident_df, 'LYSO_INFO': lyso_ident_df}

        if score_df.shape[0] > 0:

            img_name = output_folder + '\%s_%.4f_scan%.0f_%.4f_%s.png' % (pl_type, _usr_mz, ms2_spec_idx,
                                                                          _usr_ms2_rt, _usr_abbr_bulk)

            plot_spectra(_row_se, xic_dct, usr_ident_info_df,
                         _ms1_rt, _ms2_rt, ms1_df, ms2_df,
                         fa_indicator, lyso_indicator, fa_list_csv, save_img_as=img_name
                         )

    else:
        print('!!!!!!!!!!!! PR NOT in list !!!!!!!!!!!!')

    print('---------------------NEXT---------------------------')

print('>>> FINISHED!<<<')
