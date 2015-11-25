# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release.
# For more info please contact zhixu.ni@uni-leipzig.de

import ConfigParser
from LibLPPhunter.xic import XIC
from LibLPPhunter.msms import MSMS
import time


config = ConfigParser.ConfigParser()
config.read('config.ini')
st_time = time.clock()

infile_type = config.get('inputfile', 'filetype')
if infile_type.lower() == 'mzml':
    infile_name = config.get('inputfile', 'filename')
    print infile_name
    #
    xic_spec = XIC(infile_name)
    rt_lst = xic_spec.find_mz(778.560)
    print rt_lst
    msms_spec = MSMS(infile_name)
    msms_spec.get_ms2(778.560, rt_lst)

if infile_type.lower() == 'mzml.gz':
    infile_name = config.get('inputfile', 'filename')
    print infile_name

    xic_spec = XIC(infile_name)
    xic_spec.find_mz(778.560, ppm=20)

else:
    print 'No input mzML'

ed_time = time.clock() - st_time
print ed_time
