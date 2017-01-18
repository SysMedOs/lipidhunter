# -*- coding: utf-8 -*-
# Copyright 2015-2017 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

from __future__ import print_function
from __future__ import print_function
import re
import pandas as pd
import pymzml


class Extractor(object):


    def get_ms_all(self, usr_mzml, usr_th, vendor='thermo'):
        """
        _spectrum['MS:1000016'] --> Start time of the scan
        _spectrum['MS:1000744'] --> precursor for MS2
        _spectrum["MS:1000042"] --> precursor intensity
        _spectrum['MS:1000796'] --> spectrum title
        070120_CM_neg_70min_BOTH_I.1.1. File:"070120_CM_neg_70min_BOTH_I.raw", NativeID:"function=1 process=0 scan=71"
        _spectrum['MS:1000514'] --> m/z array
        _spectrum['MS:1000515'] --> intensity array

        :param usr_th:
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
        ('MS:1000511', ['value']) --> ms level

        <cvParam cvRef="MS" accession="MS:1000768" name="Thermo nativeID format" value=""/>
        <cvParam cvRef="MS" accession="MS:1000563" name="Thermo RAW format" value=""/>
        <cvParam cvRef="MS" accession="MS:1000769" name="Waters nativeID format" value=""/>
        <cvParam cvRef="MS" accession="MS:1000526" name="Waters raw format" value=""/>
        """
        waters_obo_lst = (('MS:1000016', ['value']), ('MS:1000744', ['value']), ('MS:1000042', ['value']),
                          ('MS:1000796', ['value']), ('MS:1000514', ['name']), ('MS:1000515', ['name']),
                          ('MS:1000769', ['name']), ('MS:1000526', ['name']))
        thermo_obo_lst = (('MS:1000511', ['value']), ('MS:1000768', ['name']), ('MS:1000563', ['name']))
        vendor_obo_lst = thermo_obo_lst + waters_obo_lst
        for _obo in vendor_obo_lst:
            if _obo not in pymzml.minimum.MIN_REQ:
                pymzml.minimum.MIN_REQ.append(_obo)
            else:
                pass
                # end hot patch

        _usr_spectra = pymzml.run.Reader(usr_mzml)
        _pre_ms_df = pd.DataFrame(columns=['mz', 'i', 'rt', 'scan_id'])
        _spec_title_obo = 'MS:1000796'
        _spec_level_obo = 'MS:1000511'
        _thermo_title_obo = 'MS:1000768'
        _waters_title_obo = 'MS:1000769'

        # # set default settings#
        if vendor == 'waters':
            _function_re = re.compile(r'(.*)(function=)(\d{1,2})(.*)')
            for _spectrum in _usr_spectra:

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
                            _pre_ms_df = _pre_ms_df.append(_toppeaks_df)
                    else:
                        print('NOT MS level')
                except KeyError:
                    print('Not MS', _spectrum['id'])

        elif vendor == 'thermo':
            for _spectrum in _usr_spectra:
                if _spec_level_obo in _spectrum.keys():
                    ms_level = _spectrum[_spec_level_obo]
                    print('ms_level', ms_level)
                    # _spectrum_title = _spectrum[_spec_title_obo]
                    if _spectrum[_spec_level_obo] == 1:
                        print ('MS level: %s, Scan_num: %s, Scan_time: %s;' %
                               (ms_level, _spectrum['id'], _spectrum['MS:1000016']))
                        _toppeaks_lst = _spectrum.peaks
                        _toppeaks_df = pd.DataFrame(data=_toppeaks_lst, columns=['mz', 'i'])
                        _toppeaks_df['rt'] = _spectrum['MS:1000016']
                        _toppeaks_df['scan_id'] = _spectrum['id']
                        _pre_ms_df = _pre_ms_df.append(_toppeaks_df)
                    else:
                        print('NOT MS level')

        _ms_df = _pre_ms_df[_pre_ms_df['i'] >= usr_th]
        _ms_df = _ms_df.sort_values(by='mz')
        print('_ms_df_shape', _ms_df.shape)
        return _ms_df

    def get_scan_events(self, usr_mzml, usr_ms1_abs_th, usr_ms2_abs_th, vendor = 'thermo'):

        waters_obo_lst = (('MS:1000016', ['value']), ('MS:1000744', ['value']), ('MS:1000042', ['value']),
                          ('MS:1000796', ['value']), ('MS:1000514', ['name']), ('MS:1000515', ['name']),
                          ('MS:1000769', ['name']), ('MS:1000526', ['name']))
        thermo_obo_lst = (('MS:1000511', ['value']), ('MS:1000768', ['name']), ('MS:1000563', ['name']))
        vendor_obo_lst = thermo_obo_lst + waters_obo_lst
        for _obo in vendor_obo_lst:
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
        _spec_level_obo = 'MS:1000511'
        _thermo_title_obo = 'MS:1000768'
        _waters_title_obo = 'MS:1000769'
        vendor = 'thermo'
        # for _spectrum in _usr_spectra:
        #     if _thermo_title_obo in _spectrum.keys():
        #         vendor = 'thermo'
        #         print('>>> mzML from Thermo >>>')
        #     if _waters_title_obo in _spectrum.keys():
        #         vendor = 'waters'
        #         print('>>> mzML from Waters >>>')

        # # set default settings#
        if vendor == 'waters':
            _function_re = re.compile(r'(.*)(function=)(\d{1,2})(.*)(scan=)(\d*)(.*)')
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
                            _toppeaks_df = _toppeaks_df[(_toppeaks_df['mz'] >= 400) & (_toppeaks_df['mz'] <= 1500)]
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
                        print('NOT MS level')
                except KeyError:
                    print('Not MS', _spectrum['id'])

        if vendor == 'thermo':
            for _spectrum in _usr_spectra:
                _toppeaks_lst = []
                if _spec_level_obo in _spectrum.keys():
                    rt = float(_spectrum['MS:1000016'])
                    print(rt)
                    ms_level = _spectrum[_spec_level_obo]
                    print('ms_level', ms_level)
                    # _spectrum_title = _spectrum[_spec_title_obo]
                    if ms_level == 1:
                        print ('Function: %s, Scan_num: %s, Scan_time: %s;' %
                               (ms_level, _spectrum['id'], _spectrum['MS:1000016']))
                        _toppeaks_lst = _spectrum.peaks
                        _toppeaks_df = pd.DataFrame(data=_toppeaks_lst, columns=['mz', 'i'])
                        _toppeaks_df['function'] = ms_level
                        _toppeaks_df['rt'] = _spectrum['MS:1000016']
                        _toppeaks_df['scan_id'] = _spectrum['id']
                        # _toppeaks_df = _toppeaks_df.sort_values(by='mz', ascending=False)
                        _toppeaks_df = _toppeaks_df[(_toppeaks_df['mz'] >= 400) & (_toppeaks_df['mz'] <= 1500)]
                        _toppeaks_df = _toppeaks_df[_toppeaks_df['i'] >= usr_ms1_abs_th]
                        _pre_ms_df = _pre_ms_df.append(_toppeaks_df)

                    elif ms_level == 2:
                        print ('Function: %s, Scan_num: %s, Scan_time: %s;' %
                               (ms_level, _spectrum['id'], _spectrum['MS:1000016']))
                        _toppeaks_lst = [(_spectrum['MS:1000744'], 0)]
                        _toppeaks_df = pd.DataFrame(data=_toppeaks_lst, columns=['mz', 'i'])
                        _toppeaks_df['function'] = ms_level
                        _toppeaks_df['rt'] = _spectrum['MS:1000016']
                        _toppeaks_df['scan_id'] = _spectrum['id']
                        # _toppeaks_df = _toppeaks_df.sort_values(by='mz', ascending=False)
                        # _toppeaks_df = _toppeaks_df[_toppeaks_df['i'] >= usr_ms2_abs_th]
                        _pre_ms_df = _pre_ms_df.append(_toppeaks_df)

                    else:
                        print ('NOT MS level')

        print('_ms_df_shape', _pre_ms_df.shape)
        return _pre_ms_df
