# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig          
# The software is currently  under development and is not ready to be released. 
# A suitable license will be chosen before the official release.               
# For more info please contact zhixu.ni@uni-leipzig.de


from __future__ import division

import pymzml
import pymzml.spec
import pymzml.run

import pandas as pd

class MSMS(object):

    def __init__(self, mzml):

        self.mzml_obj = pymzml.run.Reader(mzml, MS1_Precision=20e-6, MSn_Precision=200e-6)

    def get_ms2(self, mz2find, ppm=None):

        if ppm == None:
            ppm = 500
        else:
            pass

        # mz_bot = mz2find - (ppm / 1000000) * mz2find
        # mz_top = mz2find + (ppm / 1000000) * mz2find
        mz_bot = mz2find
        mz_top = mz2find + 1
        print mz_bot, mz_top

        # xic_pd = pd.DataFrame()
        # 'MS:1000511' == 'ms level'

        for spectrum in self.mzml_obj:

            try:
                if 22.0 < spectrum['MS:1000016'] < 30.0 and spectrum['MS:1000511'] == 2:
                    # print 'found1', type(spectrum['precursors']), spectrum['precursors']
                    _pr_info_dct = spectrum['precursors'][0]
                    _pr_mz = _pr_info_dct['mz']
                    # print _pr_mz
                    if mz_bot < _pr_mz < mz_top:
                        print 'found', _pr_mz, spectrum['MS:1000016']
                        # print spectrum.keys()
                        # # print spectrum['MS:1000016']
                        # print spectrum['precursors']
                        # print spectrum
                    # break

                    # Add 32-bit/64-bit encoding to mzML
                    # try:
                    #     # print spectrum['BinaryArrayOrder']
                    #     # [('encoding', '32-bit float'), ('compression', 'zlib'), ('arrayType', 'mz'),
                    #     # ('encoding', '32-bit float'), ('compression', 'zlib'), ('arrayType', 'i')]
                    #
                    #     if spectrum['BinaryArrayOrder'] == []:
                    #         spectrum['BinaryArrayOrder'] = [('encoding', '32-bit float'),
                    #                                         ('compression', 'no'), ('arrayType', 'mz'),
                    #                                         ('encoding', '32-bit float'), ('compression', 'no'),
                    #                                         ('arrayType', 'i')]
                    #         # print r"spectrum['BinaryArrayOrder'] existed"
                    #         # print spectrum['BinaryArrayOrder']
                    #     else:
                    #         # print r"spectrum['BinaryArrayOrder'] existed"
                    #         pass
                    # except KeyError:
                    #     spectrum['BinaryArrayOrder'] = [('encoding', '32-bit float'),
                    #                                     ('compression', 'no'), ('arrayType', 'mz'),
                    #                                     ('encoding', '32-bit float'), ('compression', 'no'),
                    #                                     ('arrayType', 'i')]
                    #     print 'BinaryArrayOrder Added!'

                    # print spectrum.keys()
                    # _toppeaks_lst = spectrum.highestPeaks(20)
                    # print _toppeaks_lst
                        # toppeaks_lst = []
                        # for _p in _toppeaks_lst:
                        #     p = list(_p)
                        #     p.append(spectrum['MS:1000016'])
                        #     toppeaks_lst.append(p)
                        #
                        # # print toppeaks_lst[1]
                        #
                        # _toppeaks_df = pd.DataFrame(data=toppeaks_lst, columns=['mz', 'i', 'rt'])
                        #
                        # print _toppeaks_df.head()
                        #
                        # _query_str = '( %f < mz < %f)' % (mz_bot, mz_top)
                        # # print _query_str
                        #
                        # _fund_df = _toppeaks_df.query(_query_str)
                        # _found_count = len(_fund_df['mz'].tolist())
                        # # rt_lst = [spectrum['MS:1000016']] * _found_count
                        #
                        # if _fund_df['mz'].tolist() == []:
                        #     pass
                        # else:
                        #
                        #     # _fund_df.loc[:, 'rt'] = pd.Series(rt_lst, index=_fund_df.index)
                        #     # print _fund_df.head()
                        #     xic_pd = xic_pd.append(_fund_df, ignore_index=True)

            except KeyError:
                pass

            # break


        # print xic_pd

        # xic_pd.to_csv('msms.csv')