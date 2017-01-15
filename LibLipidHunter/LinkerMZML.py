# -*- coding: utf-8 -*-
# Copyright 2015-2017 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

import re
import pandas as pd
import pymzml


def hunt_link(pl_class, usr_mzml, usr_df, params_dct):

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

    usr_ms_abs_th = params_dct['MS_THRESHOLD']
    usr_ms2_abs_th = params_dct['MS2_THRESHOLD']

    dda_top = params_dct['DDA_TOP']
    print dda_top
    # from min to sec
    rt_start = params_dct['RT_START']
    rt_end = params_dct['RT_END']
    mz_start = params_dct['MZ_START']
    mz_end = params_dct['MZ_END']

    header_basics = ['function', 'scan_id', 'mz', 'rt', 'MS1_obs_mz', 'Lib_mz', 'Abbreviation', 'Formula', 'Ion', 'Class', 'DDA_event']

    if pl_class == 'PC':
        header = header_basics + ['MZ_NL-60', 'I_NL-60', 'MZ_FRAG-168', 'I_FRAG-168', 'MZ_FRAG-224', 'I_FRAG-224']
    elif pl_class == 'SM':
        header = header_basics + ['MZ_NL-60', 'I_NL-60', 'MZ_FRAG-168', 'I_FRAG-168', 'MZ_FRAG-224', 'I_FRAG-224']
    elif pl_class == 'PE':
        header = header_basics + ['MZ_PR', 'I_PR', 'MZ_FRAG-196', 'I_FRAG-196']
    elif pl_class == 'PS':
        header = header_basics + ['MZ_NL-87', 'I_NL-87']
    elif pl_class == 'TG':
        header = header_basics
    else:
        header = header_basics

    print('header', header)

    _pre_ms_df = pd.DataFrame(columns=['function', 'scan_id', 'rt', 'mz', 'i'])
    # Waters spectra function 1 is MS level, function 2-->n-1 are MS2 spectra, function n is lock spray
    _ms2_function_lst = range(2, dda_top + 1)
    _spec_title_obo = 'MS:1000796'
    _function_re = re.compile(r'(.*)(function=)(\d{1,2})(.*)(scan=)(\d*)(.*)')

    _usr_func_lst = usr_df['function'].tolist()
    _usr_scanid_lst = usr_df['scan_id'].tolist()
    _usr_prmz_lst = usr_df['mz'].tolist()
    _usr_ident_lst = zip(_usr_func_lst, _usr_scanid_lst)

    # # set default settings#
    _out_df = pd.DataFrame()
    _ms1_df = pd.DataFrame()
    _DDA_count = 0
    for _spectrum in _usr_spectra:
        _toppeaks_lst = []

        try:
            _spectrum_title = _spectrum[_spec_title_obo]
            _function_checker = _function_re.match(_spectrum_title)
            if _function_checker:
                _function = int(_function_checker.groups()[2])
                # _scanid = int(_function_checker.groups()[5])

                if _function == 1:
                    print 'MS1 level --->>>'
                    # if _DDA_count == 0:
                    #     _DDA_count = 1
                    # else:
                    #     pass
                    _ms1_rt = _spectrum['MS:1000016']
                    if rt_start <= _ms1_rt <= rt_end:
                        _tmp_ms1_lst = _spectrum.peaks
                        _tmp_ms1_df = pd.DataFrame(data=_tmp_ms1_lst, columns=['mz', 'i'])
                        _tmp_ms1_df = _tmp_ms1_df[_tmp_ms1_df['i'] >= usr_ms_abs_th]
                        _tmp_ms1_df['DDA_event'] = _DDA_count
                        _tmp_ms1_df['rt'] = _ms1_rt
                        _ms1_df = _ms1_df.append(_tmp_ms1_df)
                        _DDA_count += 1

                if _function in _ms2_function_lst:
                    print 'MS/MS level --->>>'
                    _function = int(_function)
                    _scanid = int(_function_checker.groups()[5])
                    _prmz = float(_spectrum['MS:1000744'])
                    _prrt = float(_spectrum['MS:1000016'])
                    print (_function, _scanid, _prmz, _prrt)
                    if rt_start <= _prrt <= rt_end:
                        if (_function, _scanid) in _usr_ident_lst:
                            print 'Found scan!'
                            print ('Function: %s, Scan_num: %s, Scan_time: %f, pr_m/z: %f ;' %
                                   (_function, _scanid, _prrt, _prmz))
                            _toppeaks_lst = _spectrum.peaks
                            _toppeaks_df = pd.DataFrame(data=_toppeaks_lst, columns=['mz', 'i'])

                            _toppeaks_df = _toppeaks_df.query('i > %i' % usr_ms2_abs_th)
                            _toppeaks_df = _toppeaks_df.query('%f < mz < %f' % (mz_start, mz_end))

                            # for PC
                            if pl_class == 'PC' or pl_class == 'SM':
                                _tmp_frag168_df = _toppeaks_df.query('168 < mz < 169')
                                _tmp_frag224_df = _toppeaks_df.query('224 < mz < 225')
                                _nl_query_code = '%f < mz < %f' % (_prmz - 60.0211 - 0.5, _prmz - 60.0211 + 1)
                                _tmp_nl60_df = _toppeaks_df.query(_nl_query_code)

                                if _tmp_frag168_df.shape[0] > 0 or _tmp_frag224_df.shape[0] > 0 or _tmp_nl60_df.shape[0] > 0:
                                    print '------------------>>>>>> found key peaks'
                                    _usr_df_query_code = 'function == %i & scan_id == %i' % (_function, _scanid)
                                    # print _usr_df_query_code
                                    _tmp_usr_df = usr_df.query(_usr_df_query_code)

                                    if _tmp_nl60_df.shape[0] > 0:
                                        _tmp_nl60_df = _tmp_nl60_df.sort_values(by='i', ascending=False)
                                        _tmp_usr_df['MZ_NL-60'] = _tmp_nl60_df.iloc[0, 0]
                                        _tmp_usr_df['I_NL-60'] = _tmp_nl60_df.iloc[0, 1]
                                    if _tmp_frag168_df.shape[0] > 0:
                                        _tmp_frag168_df = _tmp_frag168_df.sort_values(by='i',
                                                                                      ascending=False)
                                        _tmp_usr_df['MZ_FRAG-168'] = _tmp_frag168_df.iloc[0, 0]
                                        _tmp_usr_df['I_FRAG-168'] = _tmp_frag168_df.iloc[0, 1]
                                    if _tmp_frag224_df.shape[0] > 0:
                                        _tmp_frag224_df = _tmp_frag224_df.sort_values(by='i',
                                                                                      ascending=False)
                                        _tmp_usr_df['MZ_FRAG-224'] = _tmp_frag224_df.iloc[0, 0]
                                        _tmp_usr_df['I_FRAG-224'] = _tmp_frag224_df.iloc[0, 1]

                                    else:
                                        pass
                                    if _tmp_usr_df.shape[0] > 0:
                                        _tmp_usr_df['rt'] = _prrt
                                        # _tmp_usr_df['prmz'] = _prmz
                                        # _tmp_usr_df['f_func'] = _function
                                        # _tmp_usr_df['f_scan'] = _scanid
                                        _tmp_usr_df['DDA_event'] = _DDA_count
                                    print('_tmp_usr_df', _tmp_usr_df.head())
                                    _out_df = _out_df.append(_tmp_usr_df)
                                    _DDA_count += 0
                                else:
                                    print 'no key peaks! --->>> skip'
                                    _DDA_count += 0

                            # for PE
                            if pl_class == 'PE':
                                _tmp_frag196_df = _toppeaks_df.query('196 < mz < 197')
                                _nl_query_code = '%f < mz < %f' % (_prmz - 0.5, _prmz + 1)
                                _tmp_pr_df = _toppeaks_df.query(_nl_query_code)

                                if _tmp_frag196_df.shape[0] > 0 or _tmp_pr_df.shape[0] > 0:
                                    print '------------------>>>>>> found key peaks'
                                    _usr_df_query_code = 'function == %i & scan_id == %i' % (_function, _scanid)
                                    # print _usr_df_query_code
                                    _tmp_usr_df = usr_df.query(_usr_df_query_code)

                                    if _tmp_pr_df.shape[0] > 0:
                                        _tmp_pr_df = _tmp_pr_df.sort_values(by='i', ascending=False)
                                        _tmp_usr_df['MZ_PR'] = _tmp_pr_df.iloc[0, 0]
                                        _tmp_usr_df['I_PR'] = _tmp_pr_df.iloc[0, 1]
                                    if _tmp_frag196_df.shape[0] > 0:
                                        _tmp_frag196_df = _tmp_frag196_df.sort_values(by='i',
                                                                                      ascending=False)
                                        _tmp_usr_df['MZ_FRAG-196'] = _tmp_frag196_df.iloc[0, 0]
                                        _tmp_usr_df['I_FRAG-196'] = _tmp_frag196_df.iloc[0, 1]

                                    else:
                                        pass
                                    if _tmp_usr_df.shape[0] > 0:
                                        _tmp_usr_df['rt'] = _prrt
                                        # _tmp_usr_df['prmz'] = _prmz
                                        # _tmp_usr_df['f_func'] = _function
                                        # _tmp_usr_df['f_scan'] = _scanid
                                        _tmp_usr_df['DDA_event'] = _DDA_count
                                    _out_df = _out_df.append(_tmp_usr_df)
                                    _DDA_count += 0
                                else:
                                    print 'no key peaks! --->>> skip'
                                    _DDA_count += 0

                            # for PS
                            if pl_class == 'PS':
                                _nl_query_code = '%f < mz < %f' % (_prmz - 87 - 0.5, _prmz - 87 + 1)
                                _tmp_nl_df = _toppeaks_df.query(_nl_query_code)

                                if _tmp_nl_df.shape[0] > 0:
                                    print '------------------>>>>>> found key peaks'
                                    _usr_df_query_code = 'function == %i & scan_id == %i' % (_function, _scanid)
                                    # print _usr_df_query_code
                                    _tmp_usr_df = usr_df.query(_usr_df_query_code)

                                    if _tmp_nl_df.shape[0] > 0:
                                        _tmp_nl_df = _tmp_nl_df.sort_values(by='i', ascending=False)
                                        _tmp_usr_df['MZ_NL-87'] = _tmp_nl_df.iloc[0, 0]
                                        _tmp_usr_df['I_NL-87'] = _tmp_nl_df.iloc[0, 1]

                                    else:
                                        pass
                                    if _tmp_usr_df.shape[0] > 0:
                                        _tmp_usr_df['rt'] = _prrt
                                        # _tmp_usr_df['prmz'] = _prmz
                                        # _tmp_usr_df['f_func'] = _function
                                        # _tmp_usr_df['f_scan'] = _scanid
                                        _tmp_usr_df['DDA_event'] = _DDA_count
                                    _out_df = _out_df.append(_tmp_usr_df)
                                    _DDA_count += 0
                                else:
                                    print 'no key peaks! --->>> skip'
                                    _DDA_count += 0

                            # for PA. PG, PIP and TG
                            else:
                                if pl_class in ['PA', 'PG', 'PI', 'PIP', 'TG']:
                                    print('for PA. PG, PIP and TG')
                                    _usr_df_query_code = 'function == %i & scan_id == %i' % (_function, _scanid)
                                    # print _usr_df_query_code
                                    _tmp_usr_df = usr_df.query(_usr_df_query_code)
                                    if _tmp_usr_df.shape[0] > 0:
                                        _tmp_usr_df['rt'] = _prrt
                                        # _tmp_usr_df['prmz'] = _prmz
                                        # _tmp_usr_df['f_func'] = _function
                                        # _tmp_usr_df['f_scan'] = _scanid
                                        _tmp_usr_df['DDA_event'] = _DDA_count
                                    _out_df = _out_df.append(_tmp_usr_df)
                                    _DDA_count += 0
                                    pass

                        else:
                            print 'not identified --->>> skip'
                            _DDA_count += 0
                else:
                    _DDA_count += 0

            else:
                print 'NOT MS level'
                _DDA_count += 0
        except KeyError:
            print 'Not MS or MS2', _spectrum['id']
    # _ms_df = _pre_ms_df.sort_values(by='mz')
    print '_out_df_shape', _out_df.shape
    print '_ms1_df_shape', _ms1_df.shape

    _out_df = _out_df.reset_index(drop=True)
    print('_out_df')
    print(_out_df.head())
    _idx_lst = _out_df.index.tolist()
    print('_idx_lst')
    print(_idx_lst)

    # for _idx in _idx_lst:
    #     print('_idx', _idx)
    #     _out_df.set_value(_idx, 'MS1_obs_mz', 0)
    #     _ms2_mz_lst = _out_df.loc[_idx, 'mz']
    #     _ms2_rt_lst = _out_df.loc[_idx, 'rt']
    #     _ms2_DDA_count_lst = _out_df.loc[_idx, 'DDA_event']
    #     print '_ms2_DDA_count_lst', _ms2_DDA_count_lst
    #
    #     try:
    #         _ms1_DDA_count = _ms2_DDA_count_lst - 1
    #         _ms1_pre_df = _ms1_df.query('DDA_event == %i' % _ms1_DDA_count)
    #         print _ms1_pre_df.shape
    #         # _ms1_pre_df = _ms1_pre_df.loc[:, ]
    #         _ms1_query = '%f < mz < %f' % (_ms2_mz_lst - 0.6, _ms2_mz_lst)
    #         _ms1_pre_df = _ms1_pre_df.query(_ms1_query)
    #         _ms1_pre_df = _ms1_pre_df.sort_values(by='i', ascending=False)
    #         print _ms1_pre_df.shape
    #         if _ms1_pre_df.shape[0] > 0:
    #             _ms1_rt = _ms1_pre_df.iloc[0, 3]
    #             _ms1_mz = _ms1_pre_df.iloc[0, 0]
    #             print _ms1_rt, _ms2_rt_lst
    #             if _ms2_rt_lst - 0.06 <= _ms1_rt < _ms2_rt_lst:
    #                 pass
    #             else:
    #                 _out_df.set_value(_idx, 'MS1_obs_mz', 0)
    #         else:
    #             print('No MS2 matches!!!')
    #         print('_ms2_DDA_count_lst is numpy.int64!')
    #
    #     except TypeError:
    #         _ms2_DDA_count_lst = _ms2_DDA_count_lst.tolist()
    #         print('_ms2_DDA_count_lst')
    #         print(_ms2_DDA_count_lst)
    #         _ms2_mz_lst = _ms2_mz_lst.tolist()
    #         _ms2_rt_lst = _ms2_rt_lst.tolist()
    #         for _ms2_DDA_count in _ms2_DDA_count_lst:
    #             _idx = _ms2_DDA_count_lst.index(_ms2_DDA_count)
    #             _ms2_mz = _ms2_mz_lst[_idx]
    #             _ms2_rt = _ms2_rt_lst[_idx]
    #             _ms1_DDA_count = _ms2_DDA_count - 1
    #             _ms1_pre_df = _ms1_df.query('DDA_event == %i' % _ms1_DDA_count)
    #             print _ms1_pre_df.shape
    #             # _ms1_pre_df = _ms1_pre_df.loc[:, ]
    #             _ms1_query = '%f < mz < %f' % (_ms2_mz - 0.6, _ms2_mz)
    #             _ms1_pre_df = _ms1_pre_df.query(_ms1_query)
    #             _ms1_pre_df = _ms1_pre_df.sort_values(by='i', ascending=False)
    #             print _ms1_pre_df.shape
    #             if _ms1_pre_df.shape[0] > 0:
    #                 _ms1_rt = _ms1_pre_df.iloc[0, 3]
    #                 _ms1_mz = _ms1_pre_df.iloc[0, 0]
    #                 print _ms1_rt, _ms2_rt
    #                 if _ms2_rt - 0.07 <= _ms1_rt < _ms2_rt:
    #                     pass
    #                 else:
    #                     _out_df.set_value(_idx, 'MS1_obs_mz', 0)
    #             else:
    #                 print('No MS2 matches!!!')
    #             print('_ms2_DDA_count_lst is pandas.Series!')
    # print(_out_df.head(5))
    _out_df['Class'] = pl_class
    final_header = []
    pre_header_lst = _out_df.columns.tolist()
    for _h in header:
        if _h in pre_header_lst:
            final_header.append(_h)
    print('final_header', final_header)
    try:
        _out_df = _out_df[final_header]
    except KeyError:
        _out_df = _out_df[header_basics]
    print(_out_df.head(5))
    return _out_df
