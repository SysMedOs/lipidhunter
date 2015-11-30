# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig          
# The software is currently  under development and is not ready to be released. 
# A suitable license will be chosen before the official release.
# For more info please contact zhixu.ni@uni-leipzig.de

from __future__ import division
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

import pymzml
import pymzml.spec
import pymzml.run

import pandas as pd

from encode_checker import add_encode_info
from LibLPPhunter.tic import TIC


class XIC(object):

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
        self.mzml = mzml

    def extract_mz(self, mz2xic_lst, rt=None, ppm=None):

        if ppm is None:
            ppm = 20
        else:
            pass

        if rt is None:
            tic_spec = TIC(self.mzml, self.encode_type)
            rt_max = tic_spec.get_maxtime()
            print rt_max
            rt = [0.1, rt_max]
        else:
            try:
                if len(rt) == 2 and rt[0] < rt[1]:
                    pass
            except:
                tic_spec = TIC(self.mzml_obj, self.encode_type)
                rt_max = tic_spec.get_maxtime()
                print rt_max
                rt = [0.1, rt_max]

        # mz_bot_lst = []
        # mz_top_lst = []
        # # mz_range_lst = []
        # for mz2xic in mz2xic_lst:
        #     mz_bot = mz2xic - (ppm / 1000000) * mz2xic
        #     mz_top = mz2xic + (ppm / 1000000) * mz2xic
        #     mz_bot_lst.append(mz_bot)
        #     mz_top_lst.append(mz_top)
        # mz_range_lst = zip(mz_bot_lst, mz_top_lst)

        # print 'Start to extract %.4f from %.2f to %.2f with %d ppm' % (mz2xic, rt[0], rt[1], ppm)

        # timeDependentIntensities = []
        '''
        ['total ion current', 'filter string', 'PY:0000000', 'id', 'MS:1000285', 'MS:1000512', 'MS:1000511',
        'MS:1000515', 'MS:1000514', 'scan start time', 'intensity array', 'defaultArrayLength', 'm/z array',
        'MS:1000128', 'ms level', 'MS:1000574', 'BinaryArrayOrder', 'profile spectrum', '32-bit float',
        'zlib compression', 'MS:1000016', 'encodedData', 'MS:1000521']

        'zlib compression': ['zlib compression', 'zlib compression']
        'MS:1000574': ['zlib compression', 'zlib compression']


        ['MS:1000511', 'MS:1000576', 'MS:1000016', 'encodedData', 'MS:1000515', None, 'MS:1000523',
        'BinaryArrayOrder', 'MS:1000514', 'defaultArrayLength', 'MS:1000127', 'PY:0000000', 'id', 'MS:1000285']

        :param mz:
        :return:
        '''
        '''
        :param mz:
        :return:
        '''
        rt_dct = {}
        xic_dct = {}
        ms_spectra_dct = {}
        for mz2xic in mz2xic_lst:
            rt_dct[mz2xic] = {}
            xic_dct[mz2xic] = pd.DataFrame()
            ms_spectra_dct[mz2xic] = {}

        # xic_df = pd.DataFrame()

        for spectrum in self.mzml_obj:
            # print spectrum.keys()
            try:
                if rt[0] < spectrum['MS:1000016'] < rt[1] and spectrum['MS:1000511'] == 1:  # 'MS:1000511' == 'ms level'

                    spectrum = add_encode_info(spectrum, self.encode_type)

                    # print spectrum.keys()
                    _toppeaks_lst = spectrum.highestPeaks(200)

                    # prepare MS spectrum

                    # todo zhixu.ni@uni-leipzig.de: rewrite this part to add rt faster.
                    toppeaks_lst = []
                    for _p in _toppeaks_lst:
                        p = list(_p)
                        p.append(spectrum['MS:1000016'])
                        toppeaks_lst.append(p)

                    # print toppeaks_lst[1]

                    _toppeaks_df = pd.DataFrame(data=toppeaks_lst, columns=['mz', 'i', 'rt'])

                    # print _toppeaks_df.head()

                    for mz2xic in mz2xic_lst:
                        mz_bot = mz2xic - (ppm / 1000000) * mz2xic
                        mz_top = mz2xic + (ppm / 1000000) * mz2xic

                        _query_str = '( %f < mz < %f)' % (mz_bot, mz_top)
                        # print _query_str

                        _found_df = _toppeaks_df.query(_query_str)
                        # _found_count = len(_found_df['mz'].tolist())
                        # rt_lst = [spectrum['MS:1000016']] * _found_count

                        if len(_found_df['mz'].tolist()) == 0:
                            pass
                        else:

                            # _found_df.loc[:, 'rt'] = pd.Series(rt_lst, index=_found_df.index)
                            # print _found_df.head()
                            print spectrum['MS:1000016'], mz2xic
                            xic_dct[mz2xic] = xic_dct[mz2xic].append(_found_df, ignore_index=True)
                            _ms_dct = ms_spectra_dct[mz2xic]
                            _mspeaks_lst = spectrum.highestPeaks(300)
                            _found_lst = (_found_df['mz'].tolist(), _found_df['i'].tolist())
                            _ms_dct[spectrum['MS:1000016']] = (_found_lst, _mspeaks_lst)


            except KeyError:
                pass

                # timeDependentIntensities.append([spectrum['MS:1000016'], I, mz])
        for mz2xic in mz2xic_lst:
            try:
                xic_dct[mz2xic] = xic_dct[mz2xic].sort_values(by='rt')
                # xic_dct[mz2xic].to_csv('xic.csv')
                rt_whole_lst = xic_dct[mz2xic]['rt'].tolist()
                rt_dct[mz2xic] = [min(rt_whole_lst), max(rt_whole_lst)]
            except:
                pass

        return rt_dct, xic_dct, ms_spectra_dct
