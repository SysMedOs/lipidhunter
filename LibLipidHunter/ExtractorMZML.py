# -*- coding: utf-8 -*-
# Copyright 2016-2017 LPP team, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of LipidHunter.
# For more info please contact:
#     LPP team oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
#     Developer Georgia Angelidou georgia.angelidou@uni-leipzig.de

from __future__ import print_function
from __future__ import print_function
import re
import pandas as pd
import pymzml


class Extractor(object):

    @staticmethod
    def get_ms_all(usr_mzml, params_dct, vendor='waters'):
        """
        _spectrum['MS:1000016'] --> Start time of the scan
        _spectrum['MS:1000744'] --> precursor for MS2
        _spectrum["MS:1000042"] --> precursor intensity
        _spectrum['MS:1000796'] --> spectrum title
        070120_CM_neg_70min_BOTH_I.1.1. File:"070120_CM_neg_70min_BOTH_I.raw", NativeID:"function=1 process=0 scan=71"
        _spectrum['MS:1000514'] --> m/z array
        _spectrum['MS:1000515'] --> intensity array

        :param params_dct:
        :param vendor:
        :param usr_mzml:
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
        _pre_ms_df = pd.DataFrame(columns=['MS2_PR_mz', 'i', 'scan_time', 'scan_number'])
        _spec_title_obo = 'MS:1000796'
        _spec_level_obo = 'MS:1000511'
        # _thermo_title_obo = 'MS:1000768'
        # _waters_title_obo = 'MS:1000769'
        # _ms2_pr_obo = 'MS:1000744'

        rt_start = params_dct['rt_start']
        rt_end = params_dct['rt_end']
        mz_start = params_dct['mz_start']
        mz_end = params_dct['mz_end']
        ms1_th = params_dct['ms1_th']

        # # set default settings#
        if vendor == 'waters':
            _function_re = re.compile(r'(.*)(function=)(\d{1,2})(.*)')
            for _spectrum in _usr_spectra:

                try:
                    _spectrum_title = _spectrum[_spec_title_obo]

                    _function_checker = _function_re.match(_spectrum_title)
                    if _function_checker:
                        _function = _function_checker.groups()[2]
                        _rt = _spectrum['MS:1000016']
                        if rt_start <= float(_rt) <= rt_end:
                            if _function == '1':
                                print ('Function: %s, Scan_num: %s, Scan_time: %s;' %
                                       (_function, _spectrum['id'], _rt))
                                _toppeaks_lst = _spectrum.peaks
                                _pre_toppeaks_df = pd.DataFrame(data=_toppeaks_lst, columns=['mz', 'i'])
                                _toppeaks_df = _pre_toppeaks_df.query('%f <= mz <= %f' % (mz_start, mz_end))
                                _toppeaks_df.loc[:, 'scan_time'] = _rt
                                _toppeaks_df.loc[:, 'scan_number'] = _spectrum['id']
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
                    _rt = _spectrum['MS:1000016']
                    if rt_start <= float(_rt) <= rt_end:
                        if _spectrum[_spec_level_obo] == 1:
                            print ('MS level: %s, Scan_num: %s, Scan_time: %s;' %
                                   (ms_level, _spectrum['id'], _rt))
                            _toppeaks_lst = _spectrum.peaks
                            _pre_toppeaks_df = pd.DataFrame(data=_toppeaks_lst, columns=['mz', 'i'])
                            _toppeaks_df = _pre_toppeaks_df.query('%f <= mz <= %f' % (mz_start, mz_end))
                            _toppeaks_df.loc[:, 'scan_time'] = _rt
                            _toppeaks_df.loc[:, 'scan_number'] = _spectrum['id']
                            _pre_ms_df = _pre_ms_df.append(_toppeaks_df)
                        else:
                            print('NOT MS level')

        _ms_df = _pre_ms_df[_pre_ms_df['i'] >= ms1_th]
        _ms_df = _ms_df.sort_values(by='mz')
        print('_ms_df_shape', _ms_df.shape)
        return _ms_df

    @staticmethod
    def get_scan_events(usr_mzml, params_dct, vendor='waters'):

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

        rt_start = params_dct['rt_start']
        rt_end = params_dct['rt_end']
        mz_start = params_dct['mz_start']
        mz_end = params_dct['mz_end']
        ms1_th = params_dct['ms1_th']
        # ms2_th = params_dct['ms2_th']
        dda_top = params_dct['dda_top']

        _usr_spectra = pymzml.run.Reader(usr_mzml)
        _pre_ms_df = pd.DataFrame(columns=['DDA_rank', 'scan_number', 'scan_time', 'MS2_PR_mz', 'i'])

        # _ms2_function_lst = ['2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13']
        _ms2_function_lst = [str(i) for i in range(2, dda_top+1)]

        _spec_title_obo = 'MS:1000796'
        _spec_level_obo = 'MS:1000511'
        # _thermo_title_obo = 'MS:1000768'
        # _waters_title_obo = 'MS:1000769'
        _ms2_pr_obo = 'MS:1000744'

        # # set default settings#
        if vendor == 'waters':
            _function_re = re.compile(r'(.*)(function=)(\d{1,2})(.*)(scan=)(\d*)(.*)')
            for _spectrum in _usr_spectra:
                # _toppeaks_lst = []
                # _rt_lst = []
                # _scan_id_lst = []
                try:
                    _spectrum_title = _spectrum[_spec_title_obo]
                    _function_checker = _function_re.match(_spectrum_title)
                    if _function_checker:
                        _function = _function_checker.groups()[2]
                        _scan = _function_checker.groups()[5]
                        _rt = _spectrum['MS:1000016']
                        if rt_start <= float(_rt) <= rt_end:

                            # if _function == '1':
                            #     print ('Function: %s, Scan_num: %s, Scan_time: %s;' %
                            #            (_function, _scan, _rt))
                            #     _toppeaks_lst = _spectrum.peaks
                            #     _toppeaks_df = pd.DataFrame(data=_toppeaks_lst, columns=['MS2_PR_mz', 'i'])
                            #     _toppeaks_df['DDA_rank'] = _function
                            #     _toppeaks_df['scan_time'] = _rt
                            #     _toppeaks_df['scan_number'] = _scan
                            #     # _toppeaks_df = _toppeaks_df.sort_values(by='mz', ascending=False)
                            #     _toppeaks_df = _toppeaks_df[(_toppeaks_df['mz'] >= mz_start) &
                            #                                 (_toppeaks_df['mz'] <= mz_end)]
                            #     _toppeaks_df = _toppeaks_df[_toppeaks_df['i'] >= ms1_th]
                            #     _pre_ms_df = _pre_ms_df.append(_toppeaks_df)

                            if _function in _ms2_function_lst:
                                if mz_start <= float(_spectrum[_ms2_pr_obo]) <= mz_end:
                                    print ('Function: %s, Scan_num: %s, Scan_time: %s, PR%s;' %
                                           (_function, _scan, _rt, _spectrum[_ms2_pr_obo]))
                                    _toppeaks_lst = [(_spectrum['MS:1000744'], 0)]
                                    _toppeaks_df = pd.DataFrame(data=_toppeaks_lst, columns=['MS2_PR_mz', 'i'])
                                    # function 1 is MS survey scan
                                    _toppeaks_df.loc[:, 'DDA_rank'] = int(_function) - 1
                                    _toppeaks_df.loc[:, 'scan_time'] = _rt
                                    _toppeaks_df.loc[:, 'scan_number'] = _scan
                                    _pre_ms_df = _pre_ms_df.append(_toppeaks_df)

                        else:
                            print('NOT in RT range')
                except KeyError:
                    print('Not MS', _spectrum['id'])

        if vendor == 'thermo':
            for _spectrum in _usr_spectra:

                dda_rank_idx = 0
                if _spec_level_obo in _spectrum.keys():
                    rt = float(_spectrum['MS:1000016'])
                    print(rt)
                    ms_level = _spectrum[_spec_level_obo]
                    print('ms_level', ms_level)
                    # _spectrum_title = _spectrum[_spec_title_obo]
                    _rt = _spectrum['MS:1000016']
                    if rt_start <= float(_rt) <= rt_end:
                        if ms_level == 1:
                            # print ('Function: %s, Scan_num: %s, Scan_time: %s;' %
                            #        (ms_level, _spectrum['id'], _rt))
                            # _toppeaks_lst = _spectrum.peaks
                            # _toppeaks_df = pd.DataFrame(data=_toppeaks_lst, columns=['MS2_PR_mz', 'i'])
                            # _toppeaks_df['DDA_rank'] = ms_level
                            # _toppeaks_df['scan_time'] = _rt
                            # _toppeaks_df['scan_number'] = _spectrum['id']
                            # # _toppeaks_df = _toppeaks_df.sort_values(by='mz', ascending=False)
                            # _toppeaks_df = _toppeaks_df[(_toppeaks_df['mz'] >= mz_start) &
                            #                             (_toppeaks_df['mz'] <= mz_end)]
                            # _toppeaks_df = _toppeaks_df[_toppeaks_df['i'] >= ms1_th]
                            # _pre_ms_df = _pre_ms_df.append(_toppeaks_df)
                            dda_rank_idx = 0

                        elif ms_level == 2:
                            dda_rank_idx += 1
                            if mz_start <= float(_spectrum[_ms2_pr_obo]) <= mz_end:
                                print ('Function: %s, Scan_num: %s, Scan_time: %s, PR%s;' %
                                       (ms_level, _spectrum['id'], _rt, _spectrum[_ms2_pr_obo]))
                                _toppeaks_lst = [(_spectrum['MS:1000744'], 0)]
                                _toppeaks_df = pd.DataFrame(data=_toppeaks_lst, columns=['MS2_PR_mz', 'i'])
                                _toppeaks_df.loc[:, 'DDA_rank'] = dda_rank_idx
                                _toppeaks_df.loc[:, 'scan_time'] = _rt
                                _toppeaks_df.loc[:, 'scan_number'] = _spectrum['id']
                                _pre_ms_df = _pre_ms_df.append(_toppeaks_df)

                    else:
                        print ('NOT in RT range')
                        dda_rank_idx = 0

        print('_ms_df_shape', _pre_ms_df.shape)
        return _pre_ms_df
