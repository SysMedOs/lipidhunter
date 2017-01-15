# -*- coding: utf-8 -*-
# Copyright 2015-2017 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

import re
import pandas as pd
import pymzml


class Extractor(object):

    def get_ms_all(self, usr_mzml, usr_th):
        """
        _spectrum['MS:1000016'] --> Start time of the scan
        _spectrum['MS:1000744'] --> precursor for MS2
        _spectrum["MS:1000042"] --> precursor intensity
        _spectrum['MS:1000796'] --> spectrum title
        070120_CM_neg_70min_BOTH_I.1.1. File:"070120_CM_neg_70min_BOTH_I.raw", NativeID:"function=1 process=0 scan=71"
        _spectrum['MS:1000514'] --> m/z array
        _spectrum['MS:1000515'] --> intensity array

        :param usr_mzml:
        :param usr_output:
        :return:
        """

        """
                hot patch for waters mzML fields
                ('MS:1000016', ['value']) --> Start time of the scan
                ('MS:1000744', ['value']) --> precursor for MS2
                ('MS:1000042', ['value']) --> precursor intensity
                ('MS:1000796', ['value']) --> spectrum title
                070120_CM_neg_70min_BOTH_I.1.1. File:"070120_CM_neg_70min_BOTH_I.raw", NativeID:"function=1 process=0 scan=71"
                ('MS:1000514', ['name']) --> m/z array
                ('MS:1000515', ['name']) --> intensity array
                """
        waters_obo_lst = (('MS:1000016', ['value']), ('MS:1000744', ['value']), ('MS:1000042', ['value']),
                          ('MS:1000796', ['value']), ('MS:1000514', ['name']), ('MS:1000515', ['name']))
        for _obo in waters_obo_lst:
            if _obo not in pymzml.minimum.MIN_REQ:
                pymzml.minimum.MIN_REQ.append(_obo)
            else:
                pass
                # end hot patch

        _usr_spectra = pymzml.run.Reader(usr_mzml)
        _pre_ms_df = pd.DataFrame(columns=['mz', 'i', 'rt', 'scan_id'])

        _spec_title_obo = 'MS:1000796'

        _function_re = re.compile(r'(.*)(function=)(\d{1,2})(.*)')

        # # set default settings#
        for _spectrum in _usr_spectra:
            _toppeaks_lst = []
            # _rt_lst = []
            # _scan_id_lst = []
            try:
                _spectrum_title = _spectrum[_spec_title_obo]
                _function_checker = _function_re.match(_spectrum_title)
                if _function_checker:
                    _function = _function_checker.groups()[2]
                    if _function == '1':
                        print ('Function: %s, Scan_num: %s, Scan_time: %s;' %
                               (_function, _spectrum['id'], _spectrum['MS:1000016']))
                        _toppeaks_lst = _spectrum.peaks
                        _toppeaks_df = pd.DataFrame(data=_toppeaks_lst, columns=['mz', 'i'])
                        _toppeaks_df['rt'] = _spectrum['MS:1000016']
                        _toppeaks_df['scan_id'] = _spectrum['id']
                        # _toppeaks_df = _toppeaks_df.sort_values(by='mz', ascending=False)
                        _pre_ms_df = _pre_ms_df.append(_toppeaks_df)
                else:
                    print 'NOT MS level'
            except KeyError:
                print 'Not MS', _spectrum['id']
        _ms_df = _pre_ms_df[_pre_ms_df['i'] >= usr_th]
        _ms_df = _ms_df.sort_values(by='mz')
        print '_ms_df_shape', _ms_df.shape
        return _ms_df

    def get_scan_events(self, usr_mzml, usr_ms1_abs_th, usr_ms2_abs_th):

        waters_obo_lst = (('MS:1000016', ['value']), ('MS:1000744', ['value']), ('MS:1000042', ['value']),
                          ('MS:1000796', ['value']), ('MS:1000514', ['name']), ('MS:1000515', ['name']))
        for _obo in waters_obo_lst:
            if _obo not in pymzml.minimum.MIN_REQ:
                pymzml.minimum.MIN_REQ.append(_obo)
            else:
                pass
                # end hot patch

        _usr_spectra = pymzml.run.Reader(usr_mzml)
        _pre_ms_df = pd.DataFrame(columns=['function', 'scan_id', 'rt', 'mz', 'i'])

        _ms2_function_lst = ['2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13']
        # _ms2_function_lst = ['2', '3', '4', '5', '6', '7']

        _spec_title_obo = 'MS:1000796'

        _function_re = re.compile(r'(.*)(function=)(\d{1,2})(.*)(scan=)(\d*)(.*)')

        # # set default settings#
        for _spectrum in _usr_spectra:
            _toppeaks_lst = []
            # _rt_lst = []
            # _scan_id_lst = []
            try:
                _spectrum_title = _spectrum[_spec_title_obo]
                _function_checker = _function_re.match(_spectrum_title)
                if _function_checker:
                    _function = _function_checker.groups()[2]
                    _scan = _function_checker.groups()[5]

                    if _function == '1':
                        print ('Function: %s, Scan_num: %s, Scan_time: %s;' %
                               (_function, _scan, _spectrum['MS:1000016']))
                        _toppeaks_lst = _spectrum.peaks
                        _toppeaks_df = pd.DataFrame(data=_toppeaks_lst, columns=['mz', 'i'])
                        _toppeaks_df['function'] = _function
                        _toppeaks_df['rt'] = _spectrum['MS:1000016']
                        _toppeaks_df['scan_id'] = _scan
                        # _toppeaks_df = _toppeaks_df.sort_values(by='mz', ascending=False)
                        _toppeaks_df = _toppeaks_df[(_toppeaks_df['mz'] >= 500) & (_toppeaks_df['mz'] <= 900)]
                        _toppeaks_df = _toppeaks_df[_toppeaks_df['i'] >= usr_ms1_abs_th]
                        _pre_ms_df = _pre_ms_df.append(_toppeaks_df)

                    if _function in _ms2_function_lst:
                        print ('Function: %s, Scan_num: %s, Scan_time: %s;' %
                               (_function, _scan, _spectrum['MS:1000016']))
                        _toppeaks_lst = [(_spectrum['MS:1000744'], 0)]
                        _toppeaks_df = pd.DataFrame(data=_toppeaks_lst, columns=['mz', 'i'])
                        _toppeaks_df['function'] = _function
                        _toppeaks_df['rt'] = _spectrum['MS:1000016']
                        _toppeaks_df['scan_id'] = _scan
                        # _toppeaks_df = _toppeaks_df.sort_values(by='mz', ascending=False)
                        # _toppeaks_df = _toppeaks_df[_toppeaks_df['i'] >= usr_ms2_abs_th]
                        _pre_ms_df = _pre_ms_df.append(_toppeaks_df)

                else:
                    print 'NOT MS level'
            except KeyError:
                print 'Not MS', _spectrum['id']
        # _ms_df = _pre_ms_df.sort_values(by='mz')
        print '_ms_df_shape', _pre_ms_df.shape
        return _pre_ms_df
