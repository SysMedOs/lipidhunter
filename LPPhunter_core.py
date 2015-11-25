# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release.
# For more info please contact zhixu.ni@uni-leipzig.de

import ConfigParser
# import pandas as pd
import time

from LibLPPhunter.xic import XIC
from LibLPPhunter.msms import MSMS
from LibLPPhunter.plot import Spectra_Ploter

from LibLPPhunter.encode_checker import check_encode


config = ConfigParser.ConfigParser()
config.read('config.ini')
st_time = time.clock()

infile_type = config.get('inputfile', 'filetype')
ms1_precision = config.get('inputfile', 'ms1_precision')
msn_precision = config.get('inputfile', 'msn_precision')

rt_range = [24.0, 30.0]

mz2get = 778.560

if infile_type.lower() == 'mzml':
    infile_name = config.get('inputfile', 'filename')
    print infile_name
    encode_typ = check_encode(infile_name)
    print 'assume to be:', encode_typ, 'encoded'

    xic_spec = XIC(infile_name, encode_typ, ms1_precision, msn_precision)

    rt_lst, xic_df, ms_spectra_dct = xic_spec.extract_mz(mz2get, rt_range)

    print 'main peaks in range', rt_lst
    print xic_df.tail(5)
    msms_spec = MSMS(infile_name, encode_typ, ms1_precision, msn_precision)
    msms_spectra_dct = msms_spec.get_ms2(mz2get, rt_lst)

    spec_plt = Spectra_Ploter()
    spec_plt.plot_all(mz2get, xic_df, ms_spectra_dct, msms_spectra_dct)

# if infile_type.lower() == 'mzml.gz':
#     infile_name = config.get('inputfile', 'filename')
#     print infile_name
#
#     xic_spec = XIC(infile_name)
#     xic_spec.extract_mz(778.560, ppm=20)

else:
    print 'No input mzML'

ed_time = time.clock() - st_time
print ed_time
