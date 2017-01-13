# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig          
# The software is currently  under development and is not ready to be released. 
# A suitable license will be chosen before the official release.               
# For more info please contact zhixu.ni@uni-leipzig.de


from __future__ import division
import re
import pandas as pd
import pymzml
import pymzml.spec
import pymzml.run

from EncodeChecker import add_encode_info


class MSMS(object):
    """
    This class is designed to operate MS/MS related processes
    """

    def __init__(self, mzml, encode_type, ms1_precision=None, msn_precision=None):
        """

        :param mzml: the file name of mzML file
        :param encode_type: `string`  the type of encoder. Can be  '32-bit float' or '64-bit float'
        :param ms1_precision: `float` the MS level tolerance
        :param msn_precision: `float` the MSn level tolerance
        :return:
        """

        waters_obo_lst = (('MS:1000016', ['value']), ('MS:1000744', ['value']), ('MS:1000042', ['value']),
                          ('MS:1000796', ['value']), ('MS:1000514', ['name']), ('MS:1000515', ['name']))
        for _obo in waters_obo_lst:
            if _obo not in pymzml.minimum.MIN_REQ:
                pymzml.minimum.MIN_REQ.append(_obo)
            else:
                pass
                # end hot patch

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

    def get_ms2(self, usr_df, ppm=1000):
        """
        Get a list of m/z, a dictinary of Retention time (RT) and ppm value to get the peak list of give m/z

        :param usr_df:
        :param ppm: ppm values, usually 5~500, default is set to 500
        :type ppm: int
        :return msms_spectra_dct: ``dict`` A dict of m/z of precursor : MS/MS peak list
        """

        waters_obo_lst = (('MS:1000016', ['value']), ('MS:1000744', ['value']), ('MS:1000042', ['value']),
                          ('MS:1000796', ['value']), ('MS:1000514', ['name']), ('MS:1000515', ['name']))
        for _obo in waters_obo_lst:
            if _obo not in pymzml.minimum.MIN_REQ:
                pymzml.minimum.MIN_REQ.append(_obo)
            else:
                pass
                # end hot patch

        if ppm is None:
            ppm = 500
        else:
            pass

        _function_lst = usr_df['function'].tolist()
        _scan_id_lst = usr_df['scan_id'].tolist()
        _pr_mz_lst = usr_df['mz'].tolist()
        _ms2_scan_lst = zip(_function_lst, _scan_id_lst)

        _pre_ms2_df = pd.DataFrame()

        msms_spectra_dct = {}
        for _mz in _pr_mz_lst:
            msms_spectra_dct[_mz] = {}

        for _spectrum in self.mzml_obj:

            _spec_title_obo = 'MS:1000796'
            _function_re = re.compile(r'(.*)(function=)(\d{1,2})(.*)(scan=)(\d*)(.*)')

            try:
                _spectrum_title = _spectrum[_spec_title_obo]
                _function_checker = _function_re.match(_spectrum_title)
                if _function_checker:
                    _function = _function_checker.groups()[2]
                    _scan = _function_checker.groups()[5]
                    print ('_ms2_function, _ms2_scan_id', _function, _scan)
                    if (int(_function), int(_scan)) in _ms2_scan_lst:
                        _idx = _ms2_scan_lst.index((int(_function), int(_scan)))
                        _pr_mz = _pr_mz_lst[_idx]
                        print ('Function: %s, Scan_num: %s, Scan_time: %s;' %
                               (_function, _scan, _spectrum['MS:1000016']))
                        _toppeaks_lst = _spectrum.highestPeaks(1000)
                        _toppeaks_df = pd.DataFrame(data=_toppeaks_lst, columns=['mz', 'i'])
                        _toppeaks_df['function'] = _function
                        _toppeaks_df['rt'] = _spectrum['MS:1000016']
                        _toppeaks_df['scan_id'] = _scan
                        # _toppeaks_df = _toppeaks_df.sort_values(by='mz', ascending=False)
                        # _toppeaks_df = _toppeaks_df[_toppeaks_df['i'] >= usr_ms2_abs_th]
                        # _pre_ms2_df = _pre_ms2_df.append(_toppeaks_df)
                        _msms = (_pr_mz, _spectrum['MS:1000016'])

                        _msms_dct = msms_spectra_dct[_pr_mz]
                        _msms_dct[_msms] = _toppeaks_df

                        # if int(_ms2_function) == 4 and int(_ms2_scan_id) == 360:
                        #     print ('_______found____750')
                        #     print (_msms)
                        #     print (_toppeaks_df)
                else:
                    print 'NOT MS level'
            except KeyError:
                print 'Not MS', _spectrum['id']

            # try:
            #     if spectrum['MS:1000511'] == 2:
            #         for mz2get in mz2get_lst:
            #             rt_lst = rt_dct[mz2get]
            #             if len(rt_lst) == 2:
            #                 rt_bot = rt_lst[0]
            #                 rt_top = rt_lst[1]
            #             else:
            #                 # print 'no rt !'
            #                 rt_bot = 0.0
            #                 rt_top = 30.0
            #             if rt_bot < spectrum['MS:1000016'] < rt_top:
            #                 _pr_info_dct = spectrum['precursors'][0]
            #                 _pr_mz = _pr_info_dct['mz']
            #
            #                 mz_bot = mz2get
            #                 mz_top = mz2get + 1
            #                 # print mz_bot, mz_top
            #
            #                 if mz_bot < _pr_mz < mz_top:
            #                     print 'found', _pr_mz, spectrum['MS:1000016']
            #
            #                     spectrum = add_encode_info(spectrum, self.encode_type)
            #
            #                     _toppeaks_lst = spectrum.highestPeaks(300)
            #                     # print _toppeaks_lst
            #                     _msms_df = pd.DataFrame()
            #                     _toppeaks_df = pd.DataFrame(data=_toppeaks_lst, columns=['mz', 'i'])
            #                     _toppeaks_df = _toppeaks_df.sort_values(by='mz', ascending=False)
            #
            #                     _msms = (_pr_mz, spectrum['MS:1000016'])
            #
            #                     _msms_dct = msms_spectra_dct[mz2get]
            #                     _msms_dct[_msms] = _toppeaks_df
            #
            # except KeyError:
            #     pass

        # for _msms in msms_spectra_dct:
        #     print(_msms)
        #     print(msms_spectra_dct[_msms])
        #     break

        return msms_spectra_dct
