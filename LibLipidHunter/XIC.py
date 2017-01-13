# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig          
# The software is currently  under development and is not ready to be released. 
# A suitable license will be chosen before the official release.
# For more info please contact zhixu.ni@uni-leipzig.de

from __future__ import print_function
from __future__ import division
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

import re

import pymzml
import pymzml.spec
import pymzml.run

import pandas as pd

from EncodeChecker import add_encode_info
from LibLipidHunter.TIC import TIC


class XIC(object):

    def __init__(self, mzml, encode_type, ms1_precision=None, msn_precision=None):

        waters_obo_lst = (('MS:1000016', ['value']), ('MS:1000744', ['value']), ('MS:1000042', ['value']),
                          ('MS:1000796', ['value']), ('MS:1000514', ['name']), ('MS:1000515', ['name']))
        for _obo in waters_obo_lst:
            if _obo not in pymzml.minimum.MIN_REQ:
                pymzml.minimum.MIN_REQ.append(_obo)
            else:
                pass

        if 0 < ms1_precision < 1:
            pass
        else:
            ms1_precision = 20e-6
        if 0 < msn_precision < 1:
            pass
        else:
            msn_precision = 200e-6
        self.mzml_obj = pymzml.run.Reader(mzml, MS1_Precision=ms1_precision, MSn_Precision=msn_precision)
        # self.encode_type = encode_type
        self.mzml = mzml

    def extract_mz(self, mz2xic_df, rt=None, ppm=None):

        waters_obo_lst = (('MS:1000016', ['value']), ('MS:1000744', ['value']), ('MS:1000042', ['value']),
                          ('MS:1000796', ['value']), ('MS:1000514', ['name']), ('MS:1000515', ['name']))
        for _obo in waters_obo_lst:
            if _obo not in pymzml.minimum.MIN_REQ:
                pymzml.minimum.MIN_REQ.append(_obo)
            else:
                pass

        if ppm is None:
            ppm = 20
        else:
            pass

        if rt is None:
            # tic_spec = TIC(self.mzml, self.encode_type)
            tic_spec = TIC(self.mzml)
            rt_max = tic_spec.get_maxtime()
            print (rt_max)
            rt = [0.1, rt_max]
        else:
            try:
                if len(rt) == 2 and rt[0] < rt[1]:
                    pass
            except IndexError:
                tic_spec = TIC(self.mzml_obj)
                # tic_spec = TIC(self.mzml_obj, self.encode_type)
                rt_max = tic_spec.get_maxtime()
                print (rt_max)
                rt = [0.1, rt_max]
        rt_l = rt[0]
        rt_r = rt[1]
        print(rt_l, rt_r)
        # mz_bot_lst = []
        # mz_top_lst = []
        # # mz_range_lst = []
        # for mz2xic in mz2xic_lst:
        #     mz_bot = mz2xic - (ppm / 1000000) * mz2xic
        #     mz_top = mz2xic + (ppm / 1000000) * mz2xic
        #     mz_bot_lst.append(mz_bot)
        #     mz_top_lst.append(mz_top)
        # mz_range_lst = zip(mz_bot_lst, mz_top_lst)

        # print 'Start to extract_mzml %.4f from %.2f to %.2f with %d ppm' % (mz2xic, rt[0], rt[1], ppm)

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
        for _i, _r in mz2xic_df.iterrows():
            mz2xic = _r['MS1_obs_mz']
            rt_dct[mz2xic] = {}
            xic_dct[mz2xic] = pd.DataFrame()
            ms_spectra_dct[mz2xic] = {}

        # xic_df = pd.DataFrame()
        spec_title_obo = 'MS:1000796'
        function_re = re.compile(r'(.*)(function=)(\d{1,2})(.*)(scan=)(\d*)(.*)')

        for spectrum in self.mzml_obj:
            spec_keys = spectrum.keys()
            # print(spec_keys)
            # spectrum = add_encode_info(spectrum, self.encode_type)
            if 'MS:1000511' in spec_keys and 'MS:1000016' in spec_keys and spec_title_obo in spec_keys:
                spectrum_title = spectrum[spec_title_obo]
                function_checker = function_re.match(spectrum_title)
                if function_checker:
                    function = int(function_checker.groups()[2])
                    scanid = int(function_checker.groups()[5])

                    scan_rt = spectrum['MS:1000016']
                    # 'MS:1000511' == 'ms level'
                    # if rt[0] < spectrum['MS:1000016'] < rt[1] and spectrum['MS:1000511'] == 1:
                    if rt_l < scan_rt < rt_r:
                        if spectrum['MS:1000511'] == 1:
                            print('Read function:', spectrum['MS:1000511'], ' @RT ', scan_rt)

                            # print spectrum.keys()
                            # _toppeaks_lst = spectrum.Peaks(500)

                            # prepare MS spectrum

                            # todo zhixu.ni@uni-leipzig.de: rewrite this part to add rt faster.
                            # toppeaks_lst = []
                            # for _p in _toppeaks_lst:
                            #     p = list(_p)
                            #     p.append(spectrum['MS:1000016'])
                            #     toppeaks_lst.append(p)
                            #
                            # # print toppeaks_lst[1]
                            #
                            # _toppeaks_df = pd.DataFrame(data=toppeaks_lst, columns=['mz', 'i', 'rt'])

                            # _toppeaks_df = pd.DataFrame(data=_toppeaks_lst, columns=['mz', 'i'])
                            # _toppeaks_df['rt'] = scan_rt
                            # print('_toppeaks_df.head()')
                            # print(_toppeaks_df.head())

                            _tmp_ms1_lst = spectrum.peaks
                            _tmp_ms1_df = pd.DataFrame(data=_tmp_ms1_lst, columns=['mz', 'i'])

                            for _i, _r in mz2xic_df.iterrows():
                                mz2xic = _r['MS1_obs_mz']
                                rt = _r['rt']

                                if rt_l <= scan_rt <= rt_r:

                                    # mz_bot = mz2xic - (ppm / 1000000) * mz2xic
                                    # mz_top = mz2xic + (ppm / 1000000) * mz2xic
                                    mz_bot = mz2xic - 0.0001
                                    mz_top = mz2xic + 0.0001

                                    _query_str = '(%f <= mz <= %f)' % (mz_bot, mz_top)

                                    _found_df = _tmp_ms1_df.query(_query_str)
                                    _found_df['rt'] = scan_rt
                                    # _found_df = _toppeaks_df.query(_query_str)
                                    # _found_count = len(_found_df['mz'].tolist())
                                    # rt_lst = [spectrum['MS:1000016']] * _found_count

                                    if _found_df.shape[0] == 0:
                                        pass
                                    else:
                                        xic_dct[mz2xic] = xic_dct[mz2xic].append(_found_df, ignore_index=True)
                                        _ms2_pr_th_chker_df = _found_df.query('1000 <= i')
                                        _ms_dct = ms_spectra_dct[mz2xic]

                                        if _ms2_pr_th_chker_df.shape[0] > 0:
                                            _tmp_ms1_df = _tmp_ms1_df.sort_values(by='i', ascending=False)
                                            _tmp_ms1_df = _tmp_ms1_df.head(1000)
                                            # print(_tmp_ms1_df.head())
                                            _ms_dct[scan_rt] = (_found_df, _tmp_ms1_df, function, scanid)
                                        # _found_df.loc[:, 'rt'] = pd.Series(rt_lst, index=_found_df.index)
                                        # print _found_df.head()
                                        # print(scan_rt, mz2xic)
                                        # print(_query_str)
                                        # xic_dct[mz2xic] = xic_dct[mz2xic].append(_found_df, ignore_index=True)
                                        # _ms_dct = ms_spectra_dct[mz2xic]
                                        # _mspeaks_lst = spectrum.highestPeaks(500)
                                        # _found_lst = (_found_df['mz'].tolist(), _found_df['i'].tolist())
                                        # print('_found_lst')
                                        # print(_found_lst)

                                        # _mspeaks_lst = spectrum.peaks
                                        # # Trigger MS/MS at MS precursor > 1000
                                        # _found_lst = []
                                        # for _pre_found_mz in _mspeaks_lst:
                                        #     if _pre_found_mz[1] > 1000:
                                        #         _found_lst.append(_pre_found_mz)
                                        # if rt - 1 <= spectrum['MS:1000016'] <= rt:
                                        # _ms_dct[scan_rt] = (_found_lst, _mspeaks_lst)
            else:
                print('not MS1')
        # for _i, _r in mz2xic_df.iterrows():
        #     mz2xic = _r['MS1_obs_mz']
        #     try:
        #         xic_dct[mz2xic] = xic_dct[mz2xic].sort_values(by='rt')
        #         # xic_dct[mz2xic].to_csv('xic.csv')
        #         rt_whole_lst = xic_dct[mz2xic]['rt'].tolist()
        #         rt_dct[mz2xic] = [min(rt_whole_lst), max(rt_whole_lst)]
        #     except:
        #         pass

        return rt_dct, xic_dct, ms_spectra_dct
