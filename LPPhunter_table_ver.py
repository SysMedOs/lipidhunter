# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release.
# For more info please contact zhixu.ni@uni-leipzig.de

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
infile_name = r'D:\LPPhunter\spectra_mzML\070120_CM_neg_70min_SIN_II.mzML'
usr_xlsx = r'D:\LPPhunter\Marked_PC_70min_mgfHunter_3.xlsx'
output_folder = 'D:\LPPhunter\images'

ms1_precision = 100e-6
msn_precision = 200e-6
rt_range = [20, 30]
usr_df = pd.read_excel(usr_xlsx)

print usr_df.head()

mz1get_lst = usr_df['MS1_obs_mz'].tolist()
mz1get_lst.sort()
mz2get_lst = usr_df['mz'].tolist()
mz2get_lst.sort()
print (mz2get_lst[0:10])

print len(mz2get_lst), ' different m/z need to find===>'


# for mz2get in mz2get_lst:

print infile_name
encode_typ = check_encode(infile_name)
print 'assume to be:', encode_typ, 'encoded'

xic_spec = XIC(infile_name, encode_typ, ms1_precision, msn_precision)

rt_dct, xic_dct, ms_spectra_dct = xic_spec.extract_mz(mz1get_lst, rt_range)


msms_spec = MSMS(infile_name, encode_typ, ms1_precision, msn_precision)
msms_spectra_dct = msms_spec.get_ms2(usr_df)

# print ('msms_spectra_dct')
# print (msms_spectra_dct)
# print ('msms_spectra_dct')
# print (msms_spectra_dct)
spec_plt = Spectra_Ploter()
spec_plt.plot_all(mz1get_lst, mz2get_lst, usr_df, xic_dct, ms_spectra_dct, msms_spectra_dct, path=output_folder)

ed_time = time.clock() - st_time
print ed_time
