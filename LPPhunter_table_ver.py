# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release.
# For more info please contact zhixu.ni@uni-leipzig.de

from __future__ import print_function
import ConfigParser
import time

import pandas as pd

from rdkit import Chem

from LibLPPhunter.XIC import XIC
from LibLPPhunter.MSMS import MSMS
from LibLPPhunter.Plot_table_ver import Spectra_Ploter

from LibLPPhunter.EncodeChecker import check_encode
from LibLPPhunter.ExactMassCalc import Elem2Mass


st_time = time.clock()
print('Start --->')
Pl_type = 'PS'
# sp_name = '70min_SIN_II'
# infile_name = r'D:\LPPhunter\spectra_mzML\070120_CM_neg_%s.mzML' % sp_name
# usr_xlsx = r'D:\LPPhunter\inputfiles\%s\Marked_%s_70min_mgfHunter_3.xlsx' % (sp_name, pl_type)
# output_folder = r'D:\LPPhunter\images\%s\%s_3' % (sp_name, pl_type)
# fa_list_csv = r'D:\LPPhunter\FA_list.csv'
infile_name = r'D:\project_mzML\CM_DDA_neg_mzML\070120_CM_neg_70min_SIN_II.mzML'
usr_xlsx = r'D:\project_mzML\CM_DDA_neg_mzML\extractor_output\PC\PC_70min_SIN_II.xlsx'
output_folder = r'D:\project_mzML\CM_DDA_neg_mzML\images\%s\70min_SIN_II' % Pl_type
fa_list_csv = r'D:\LPPhunter\FA_list.csv'

ms1_precision = 0.5
msn_precision = 500e-6
rt_range = [20, 31]
usr_df = pd.read_excel(usr_xlsx)

print(usr_df.head())

# mz1get_lst = usr_df['MS1_obs_mz'].tolist()
# mz1get_lst.sort()

mz1_get_df = usr_df[['MS1_obs_mz', 'rt']]
mz1_get_df = mz1_get_df.sort_values(by=['MS1_obs_mz', 'rt'])
mz1_get_df = mz1_get_df.query('%f<= rt <= %f' % (rt_range[0], rt_range[1]))
usr_df = usr_df.sort_values(by=['mz'])
mz1_get_lst = usr_df['MS1_obs_mz'].tolist()
# print(mz1_get_lst)
mz2get_lst = usr_df['mz'].tolist()
# mz2get_lst.sort()
# print (mz2get_lst[0:10])

print(len(mz2get_lst), ' different m/z need to find===>')


# for mz2get in mz2get_lst:

print(infile_name)
encode_typ = check_encode(infile_name)
print('assume to be:', encode_typ, 'encoded')

xic_spec = XIC(infile_name, encode_typ, ms1_precision, msn_precision)

rt_dct, xic_dct, ms_spectra_dct = xic_spec.extract_mz(mz1_get_df, rt_range)
# _test_df = ms_spectra_dct[734.529541015625][25.0342999][1]
# _test_df.to_csv('_test_df')
# print(ms_spectra_dct[734.529541015625][25.0342999][0])


msms_spec = MSMS(infile_name, encode_typ, ms1_precision, msn_precision)
msms_spectra_dct = msms_spec.get_ms2(usr_df, ppm=1000)

# print ('msms_spectra_dct')
# print (msms_spectra_dct)
spec_plt = Spectra_Ploter()

spec_plt.plot_all(mz1_get_lst, mz2get_lst, usr_df, xic_dct, ms_spectra_dct,
                  msms_spectra_dct, fa_list_csv, Pl_type, path=output_folder)

ed_time = time.clock() - st_time
print(ed_time)
