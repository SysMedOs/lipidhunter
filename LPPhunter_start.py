# -*- coding: utf-8 -*-
# Copyright 2015-2016 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de


import time

import pandas as pd

from rdkit import Chem

from LibLPPhunter.XIC import XIC
from LibLPPhunter.MSMS import MSMS
from LibLPPhunter.Plot_table_ver import Spectra_Ploter

from LibLPPhunter.EncodeChecker import check_encode
from LibLPPhunter.ExactMassCalc import Elem2Mass

from LibLPPhunter.SpectraExtractor import extract_mzml
from LibLPPhunter.SpectraExtractor import get_spectra
from LibLPPhunter.SpectraExtractor import get_xic
from LibLPPhunter.PanelPloter import plot_spectra

from LibLPPhunter.FAindicator import FAindicator
from LibLPPhunter.FAindicator import Lyso_indicator

st_time = time.clock()

pl_type = 'PE'
usr_mzml = r'D:\project_mzML\CM_DDA_neg_mzML\070120_CM_neg_30min_SIN_II.mzML'
usr_xlsx = r'D:\project_mzML\CM_DDA_neg_mzML\extractor_output\%s\%s_30min_SIN_II.xlsx' % (pl_type, pl_type)
output_folder = r'D:\project_mzML\CM_DDA_neg_mzML\images\%s\30min_SIN_II' % pl_type
fa_list_csv = r'D:\LPPhunter\FA_list.csv'

usr_rt_range = [24, 25]
usr_pr_mz_range = [600, 1000]
usr_dda_top = 12
usr_ms1_precision = 50e-6
usr_msn_precision = 500e-6

print('=== ==> Start to process')
print('=== ==> --> Phospholipid class: %s' % pl_type)

usr_df = pd.read_excel(usr_xlsx)
usr_df = usr_df.round({'mz': 6})
usr_df = usr_df.query('%f<= rt <= %f' % (usr_rt_range[0], usr_rt_range[1]))
usr_df = usr_df.query('%f<= mz <= %f' % (usr_pr_mz_range[0], usr_pr_mz_range[1]))
usr_df = usr_df.sort_values(by=['rt'])
mz2get_lst = usr_df['mz'].tolist()

print('=== ==> --> Total precursor number: %i' % usr_df.shape[0])

fa_indicator = FAindicator(fa_list_csv)
lyso_indicator = Lyso_indicator(fa_list_csv)

# extract all spectra from mzML to pandas DataFrame
usr_scan_info_df, usr_spectra_pl = extract_mzml(usr_mzml, usr_rt_range, dda_top=usr_dda_top,
                                                ms1_precision=usr_ms1_precision, msn_precision=usr_msn_precision
                                                )

# get spectra of one ABBR and plot
for _i, _row_se in usr_df.iterrows():
    print(_row_se)
    _usr_mz = _row_se['mz']
    _usr_ms1_obs_mz = _row_se['MS1_obs_mz']
    _usr_ms2_rt = _row_se['rt']
    _usr_abbr_bulk = _row_se['Abbreviation']
    _usr_formula = _row_se['Formula']
    _usr_ms2_function = _row_se['function']
    _usr_ms2_scan_id = _row_se['scan_id']
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

        _xic_df = get_xic(_usr_ms1_obs_mz, usr_mzml, usr_rt_range,
                          ppm=50, ms1_precision=50e-6, msn_precision=500e-6)

        _usr_abbr_bulk = _usr_abbr_bulk.replace('(', '[')
        _usr_abbr_bulk = _usr_abbr_bulk.replace(')', ']')
        _usr_abbr_bulk = _usr_abbr_bulk.replace(':', '-')
        _usr_abbr_bulk = _usr_abbr_bulk.replace('\\', '_')
        _usr_abbr_bulk = _usr_abbr_bulk.replace('/', '_')

        print ('_usr_abbr_bulk', _usr_abbr_bulk)
        img_name = output_folder + '\%s_%.4f_scan%.0f_%.4f_%s.png' % (pl_type, _usr_mz, ms2_spec_idx,
                                                                      _usr_ms2_rt, _usr_abbr_bulk)

        plot_spectra(_row_se, _xic_df,
                     _ms1_rt, _ms2_rt, ms1_df, ms2_df,
                     fa_indicator, lyso_indicator, fa_list_csv, save_img_as=img_name
                     )

    else:
        print('!!!!!!!!!!!! PR NOT in list !!!!!!!!!!!!')

    print('---------------------NEXT---------------------------')

print('>>> FINISHED!<<<')
