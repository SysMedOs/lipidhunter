# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig          
# The software is currently  under development and is not ready to be released. 
# A suitable license will be chosen before the official release.               
# For more info please contact zhixu.ni@uni-leipzig.de

from __future__ import division
import pymzml

from encode_checker import add_encode_info


class TIC(object):
    def __init__(self, mzml, encode_type, ms1_precision=None, msn_precision=None):
        if 0 < ms1_precision < 1:
            pass
        else:
            ms1_precision = 20e-6
        if 0 < msn_precision < 1:
            pass
        else:
            msn_precision = 200e-6
        self.mzml_obj = pymzml.run.Reader(mzml, MS1_Precision=ms1_precision, MSn_Precision=msn_precision)
        self.encode_type = encode_type

    def get_maxtime(self):

        rt_max = 0.0

        for spectrum in self.mzml_obj:
            if spectrum['id'] != 'TIC':
                if spectrum['MS:1000511'] == 1:
                    if spectrum['MS:1000016'] > rt_max:
                        rt_max = spectrum['MS:1000016']

        return rt_max

