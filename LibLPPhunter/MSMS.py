# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig          
# The software is currently  under development and is not ready to be released. 
# A suitable license will be chosen before the official release.               
# For more info please contact zhixu.ni@uni-leipzig.de


from __future__ import division
import pandas as pd
import pymzml
import pymzml.spec
import pymzml.run

from EncodeChecker import add_encode_info

class MSMS(object):

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

    def get_ms2(self, mz2get_lst, rt_dct, ppm=None):

        # if len(rt_lst) == 2:
        #     rt_bot = rt_lst[0]
        #     rt_top = rt_lst[1]
        #
        # else:
        #     print 'no rt'
        #     rt_bot = 0.0
        #     rt_top = 30.0

        if ppm == None:
            ppm = 500
        else:
            pass

        msms_spectra_dct = {}
        for mz2get in mz2get_lst:
            msms_spectra_dct[mz2get] = {}

        for spectrum in self.mzml_obj:

            try:
                if spectrum['MS:1000511'] == 2:
                    for mz2get in mz2get_lst:
                        rt_lst = rt_dct[mz2get]
                        if len(rt_lst) == 2:
                            rt_bot = rt_lst[0]
                            rt_top = rt_lst[1]
                        else:
                            print 'no rt !'
                            rt_bot = 0.0
                            rt_top = 30.0
                        if rt_bot < spectrum['MS:1000016'] < rt_top:
                            _pr_info_dct = spectrum['precursors'][0]
                            _pr_mz = _pr_info_dct['mz']

                            mz_bot = mz2get
                            mz_top = mz2get + 1
                            print mz_bot, mz_top

                            if mz_bot < _pr_mz < mz_top:
                                print 'found', _pr_mz, spectrum['MS:1000016']

                                spectrum = add_encode_info(spectrum, self.encode_type)

                                _toppeaks_lst = spectrum.highestPeaks(300)
                                # print _toppeaks_lst
                                _msms_df = pd.DataFrame()
                                _toppeaks_df = pd.DataFrame(data=_toppeaks_lst, columns=['mz', 'i'])
                                _toppeaks_df = _toppeaks_df.sort_values(by='mz', ascending=False)

                                _msms = (_pr_mz, spectrum['MS:1000016'])

                                _msms_dct = msms_spectra_dct[mz2get]
                                _msms_dct[_msms] = _toppeaks_df

            except KeyError:
                pass

        # for _msms in msms_spectra_dct:
        #     print(_msms)
        #     print(msms_spectra_dct[_msms])
        #     break

        return msms_spectra_dct
