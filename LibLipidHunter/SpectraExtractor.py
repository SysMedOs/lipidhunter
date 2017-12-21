# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2017  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
# LipidHunter is Dual-licensed
#     For academic and non-commercial use: `GPLv2 License` Please read more information by the following link:
#         [The GNU General Public License version 2] (https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
#     For commercial use:
#         please contact the SysMedOs_team by email.
# Please cite our publication in an appropriate form.
# Ni, Zhixu, Georgia Angelidou, Mike Lange, Ralf Hoffmann, and Maria Fedorova.
# "LipidHunter identifies phospholipids by high-throughput processing of LC-MS and shotgun lipidomics datasets."
# Analytical Chemistry (2017).
# DOI: 10.1021/acs.analchem.7b01126
#
# For more info please contact:
#     SysMedOs_team: oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
#     Developer Georgia Angelidou georgia.angelidou@uni-leipzig.de

from __future__ import division
from __future__ import print_function
import re

import pandas as pd
import pymzml

from ParallelFunc import ppm_window_para


def extract_mzml(mzml, rt_range, dda_top=6, ms1_threshold=1000, ms2_threshold=10,
                 ms1_precision=50e-6, ms2_precision=500e-6, vendor='waters', ms1_max=0):

    """
    Extract mzML to a scan info DataFrame and a pandas panel for spectra DataFrame of mz and i

    :param ms1_max: Max of MS1 intensity, required for better identification for low intensity signals
    :param vendor: 'waters' or 'thermo'
    :param ms2_threshold:
    :param ms1_threshold:
    :param mzml: The file path of mzML file
    :type mzml: str
    :param rt_range: A List of RT. e.g. [15, 30] for 15 to 30 min
    :type rt_range: list
    :param dda_top: DDA settings e.g. DDA TOP 6
    :type dda_top: int
    :param ms1_precision: e.g. 50e-6 for 50 ppm
    :type ms1_precision: float
    :param ms2_precision: e.g. 500e-6 for 500 ppm
    :type ms2_precision: float
    :returns: scan_info_df, spec_pl
    :rtype: pandas.DataFrame, pandas.Panel

    """

    waters_obo_lst = [('MS:1000016', ['value']), ('MS:1000744', ['value']), ('MS:1000042', ['value']),
                      ('MS:1000796', ['value']), ('MS:1000514', ['name']), ('MS:1000515', ['name']),
                      ('MS:1000769', ['name']), ('MS:1000526', ['name']),
                      ('MS:1002476', ['name']), ('MS:1000514', ['name']), ('MS:1000515', ['name'])]
    thermo_obo_lst = [('MS:1000511', ['value']), ('MS:1000768', ['name']), ('MS:1000563', ['name'])]
    vendor_obo_lst = thermo_obo_lst + waters_obo_lst
    for _obo in vendor_obo_lst:
        if _obo not in pymzml.minimum.MIN_REQ:
            pymzml.minimum.MIN_REQ.append(_obo)
        else:
            pass
            # end hot patch

    rt_start = rt_range[0]
    rt_end = rt_range[1]

    print('==> Start to process file: %s' % mzml)
    print('=== ==> RT: %.2f -> %.2f with DDA Top % i' % (rt_start, rt_end, dda_top))

    spec_title_obo = 'MS:1000796'
    scan_rt_obo = 'MS:1000016'
    scan_pr_mz_obo = 'MS:1000744'
    _spec_level_obo = 'MS:1000511'

    ims_obo = 'MS:1002476'
    mz_array_obo = 'MS:1000514'
    i_array_obo = 'MS:1000515'

    spec_obj = pymzml.run.Reader(mzml, MS1_Precision=ms1_precision, MSn_Precision=ms2_precision)

    spec_idx = 0
    dda_event_idx = 0
    spec_idx_lst = []
    dda_event_lst = []
    rt_lst = []
    dda_rank_lst = []
    scan_id_lst = []
    pr_mz_lst = []

    scan_info_dct = {'spec_index': spec_idx_lst, 'scan_time': rt_lst, 'dda_event_idx': dda_event_lst,
                     'DDA_rank': dda_rank_lst, 'scan_number': scan_id_lst, 'MS2_PR_mz': pr_mz_lst}

    spec_dct = {}

    ms2_function_range_lst = range(2, dda_top + 1)
    function_range_lst = range(1, dda_top + 1)

    ms1_xic_df = pd.DataFrame()

    print('Instrument vendor: %s' % vendor)

    dda_add = 0

    if vendor == 'waters':
        scan_info_re = re.compile(r'(.*)(function=)(\d{1,2})(.*)(scan=)(\d*)(.*)')
        for _spectrum in spec_obj:
            pr_mz = 0

            not_empty_spec = 0
            if ims_obo in _spectrum.keys():
                try:
                    if len(_spectrum.peaks) > 2:
                        not_empty_spec = 1
                    else:
                        not_empty_spec = 0
                except:
                    not_empty_spec = 0
            else:
                not_empty_spec = 1

            if not_empty_spec == 1 and spec_title_obo in _spectrum.keys() and scan_rt_obo in _spectrum.keys():
                _spectrum_title = _spectrum[spec_title_obo]
                _scan_rt = float(_spectrum[scan_rt_obo])
                scan_info_checker = scan_info_re.match(_spectrum_title)

                if rt_start <= _scan_rt <= rt_end and scan_info_checker:
                    _function = int(scan_info_checker.groups()[2])
                    _scan_id = int(scan_info_checker.groups()[5])
                    if _function in function_range_lst:

                        if _function == 1:

                            _tmp_spec_df = pd.DataFrame(data=_spectrum.peaks, columns=['mz', 'i'])
                            # if _tmp_spec_df.shape[0] > 0:
                            # dda_event_idx += 1
                            dda_add = 1

                            # _tmp_spec_df = _tmp_spec_df.sort_values(by='i', ascending=False).head(1000)
                            if ms1_max > ms1_threshold:
                                _tmp_spec_df = _tmp_spec_df.query('%f <= i <= %f' %
                                                                  ((ms1_threshold * 0.1), ms1_max))
                            else:
                                _tmp_spec_df = _tmp_spec_df.query('%f <= i' % (ms1_threshold * 0.1))
                            if _tmp_spec_df.shape[0] > 0:
                                _tmp_spec_df = _tmp_spec_df.sort_values(by='i', ascending=False)
                                _tmp_spec_df = _tmp_spec_df.reset_index(drop=True)
                                spec_dct[spec_idx] = _tmp_spec_df
                                # _tmp_spec_df.loc[:, 'rt'] = _scan_rt
                                for _idx, _r in _tmp_spec_df.iterrows():
                                    _tmp_spec_df.set_value(_idx, 'rt', _scan_rt)
                                ms1_xic_df = ms1_xic_df.append(_tmp_spec_df)
                            # else:
                            #     print('empty_MS1_spectrum --> index = ', spec_idx)
                            del _tmp_spec_df

                            spec_idx_lst.append(spec_idx)
                            dda_event_lst.append(dda_event_idx)
                            rt_lst.append(_scan_rt)
                            dda_rank_lst.append(_function - 1)  # function 1 in Waters file is MS level
                            scan_id_lst.append(_scan_id)
                            pr_mz_lst.append(pr_mz)
                            print('MS1_spectrum --> index = ', spec_idx,  '; DDA_events = ', dda_event_idx)

                        if _function in ms2_function_range_lst:

                            _tmp_spec_df = pd.DataFrame(data=_spectrum.peaks, columns=['mz', 'i'])
                            pr_mz = _spectrum[scan_pr_mz_obo]
                            _ms2_temp_spec_df = _tmp_spec_df.query('i >= %f' % ms2_threshold)
                            if _ms2_temp_spec_df.shape[0] > 0:
                                spec_dct[spec_idx] = _ms2_temp_spec_df
                                del _ms2_temp_spec_df
                            # else:
                            #     print('empty_MS2_spectrum --> index = ', spec_idx)

                            del _tmp_spec_df

                            spec_idx_lst.append(spec_idx)
                            dda_event_lst.append(dda_event_idx)
                            rt_lst.append(_scan_rt)
                            dda_rank_lst.append(_function - 1)  # function 1 in Waters file is MS level
                            scan_id_lst.append(_scan_id)
                            pr_mz_lst.append(pr_mz)

                            if dda_add == 1:
                                dda_event_idx += 1
                                dda_add = 0
                            else:
                                pass
                            print('MS2_spectrum --> index = ', spec_idx, 'pr_mz:', pr_mz, 'scan_time:', _scan_rt,
                                  '; DDA_events = ', dda_event_idx)

            spec_idx += 1

    elif vendor == 'thermo':
        dda_rank_idx = 0
        for _spectrum in spec_obj:
            pr_mz = 0
            if spec_title_obo in _spectrum.keys() and scan_rt_obo in _spectrum.keys():
                # _spectrum_title = _spectrum[spec_title_obo]
                _scan_rt = float(_spectrum[scan_rt_obo])
                ms_level = _spectrum[_spec_level_obo]
                if rt_start <= _scan_rt <= rt_end:
                    _scan_id = _spectrum['id']
                    if ms_level in function_range_lst:
                        _tmp_spec_df = pd.DataFrame(data=_spectrum.peaks, columns=['mz', 'i'])
                        if ms_level == 1:
                            dda_event_idx += 1
                            dda_rank_idx = 0

                            # _tmp_spec_df = _tmp_spec_df.sort_values(by='i', ascending=False).head(1000)
                            if ms1_max > ms1_threshold:
                                _tmp_spec_df = _tmp_spec_df.query('%f <= i <= %f' % ((ms1_threshold * 0.1), ms1_max))
                            else:
                                _tmp_spec_df = _tmp_spec_df.query('%f <= i' % (ms1_threshold * 0.1))
                            if _tmp_spec_df.shape[0] > 0:
                                _tmp_spec_df = _tmp_spec_df.sort_values(by='i', ascending=False)
                                _tmp_spec_df = _tmp_spec_df.reset_index(drop=True)
                                spec_dct[spec_idx] = _tmp_spec_df
                                # _tmp_spec_df.loc[:, 'rt'] = _scan_rt
                                print('Reading MS1_survey_scan @:', _scan_rt)
                                for _idx, _r in _tmp_spec_df.iterrows():
                                    _tmp_spec_df.set_value(_idx, 'rt', _scan_rt)
                                ms1_xic_df = ms1_xic_df.append(_tmp_spec_df)
                            else:
                                print('empty_MS1_spectrum --> index = ', spec_idx)
                            del _tmp_spec_df

                        if ms_level == 2:
                            dda_rank_idx += 1
                            pr_mz = _spectrum[scan_pr_mz_obo]
                            _ms2_temp_spec_df = _tmp_spec_df.query('i >= %f' % ms2_threshold)
                            if _ms2_temp_spec_df.shape[0] > 0:
                                spec_dct[spec_idx] = _ms2_temp_spec_df
                                del _ms2_temp_spec_df
                            else:
                                print('empty_MS2_spectrum --> index = ', spec_idx)

                            del _tmp_spec_df

                        spec_idx_lst.append(spec_idx)
                        dda_event_lst.append(dda_event_idx)
                        rt_lst.append(_scan_rt)
                        dda_rank_lst.append(dda_rank_idx)
                        scan_id_lst.append(_scan_id)
                        pr_mz_lst.append(pr_mz)

                        spec_idx += 1

    scan_info_df = pd.DataFrame(data=scan_info_dct, columns=['dda_event_idx', 'spec_index', 'scan_time',
                                                             'DDA_rank', 'scan_number', 'MS2_PR_mz']
                                )

    scan_info_df = scan_info_df.sort_values(by='scan_time')
    scan_info_df = scan_info_df.round({'MS2_PR_mz': 6})

    # scan_info_df.to_excel('scan_info_df.xlsx')

    spec_pl = pd.Panel(data=spec_dct)
    print('=== ==> --> mzML extracted')

    # ms1_xic_df = ms1_xic_df.query('%f <= rt <= %f' % (rt_start, rt_end))

    return scan_info_df, spec_pl, ms1_xic_df


def get_spectra(mz, mz_lib, func_id, ms2_scan_id, ms1_obs_mz_lst,
                scan_info_df, spectra_pl, dda_top=12, ms1_precision=50e-6, vendor='waters'):

    # ms1_pr_se = pd.Series()
    ms1_df = pd.DataFrame()
    ms2_df = pd.DataFrame()
    ms1_spec_idx = 0
    ms2_spec_idx = 0
    ms1_rt = 0
    ms2_rt = 0
    ms1_mz = 0
    ms1_i = 0
    ms1_pr_ppm = 0
    function_max = dda_top + 1

    if mz in scan_info_df['MS2_PR_mz'].tolist():
        _tmp_mz_scan_info_df = scan_info_df.query('MS2_PR_mz == %.6f and DDA_rank == %f and scan_number == %f'
                                                  % (mz, func_id, ms2_scan_id)
                                                  )

        if _tmp_mz_scan_info_df.shape[0] == 1:
            ms2_spec_idx = _tmp_mz_scan_info_df.get_value(_tmp_mz_scan_info_df.index[0], 'spec_index')
            ms2_dda_idx = _tmp_mz_scan_info_df.get_value(_tmp_mz_scan_info_df.index[0], 'dda_event_idx')
            ms2_function = _tmp_mz_scan_info_df.get_value(_tmp_mz_scan_info_df.index[0], 'DDA_rank')
            ms2_scan_id = _tmp_mz_scan_info_df.get_value(_tmp_mz_scan_info_df.index[0], 'scan_number')
            ms2_rt = _tmp_mz_scan_info_df.get_value(_tmp_mz_scan_info_df.index[0], 'scan_time')

            print('%.6f @ DDA#: %.0f | Total scan id: %.0f | DDA_Rank: %.0f | Scan ID: %.0f | RT: %.4f'
                  % (mz, ms2_dda_idx, ms2_spec_idx, ms2_function, ms2_scan_id, ms2_rt)
                  )

            # get spectra_df of corresponding MS survey scan
            tmp_ms1_info_df = scan_info_df.query('dda_event_idx == %i and DDA_rank == 0' % ms2_dda_idx)
            if tmp_ms1_info_df.shape[0] > 0 and ms2_function <= function_max:
                ms1_spec_idx = tmp_ms1_info_df['spec_index'].tolist()[0]
                ms1_rt = tmp_ms1_info_df['scan_time'].tolist()[0]
                if ms1_spec_idx in spectra_pl.items:
                    ms1_df = spectra_pl[ms1_spec_idx]
                    ms1_df = ms1_df.query('i > 0')
                    ms1_df = ms1_df.sort_values(by='i', ascending=False).reset_index(drop=True)
                    ms1_delta = mz_lib * ms1_precision
                    if vendor == 'thermo':
                        ms1_pr_query = '%.7f <= mz <= %.7f' % (mz_lib - ms1_delta, mz_lib + ms1_delta)
                    else:
                        ms1_pr_query = '%.6f <= mz <= %.6f' % (mz_lib - ms1_delta, mz_lib + ms1_delta)

                    ms1_pr_df = ms1_df.query(ms1_pr_query)
                    if ms1_pr_df.shape[0] > 0:
                        ms1_pr_df.loc[:, 'mz_xic'] = ms1_pr_df['mz']
                        ms1_pr_df = ms1_pr_df.round({'mz': 6, 'mz_xic': 4})
                        ms1_pr_df = ms1_pr_df[ms1_pr_df['mz_xic'].isin(ms1_obs_mz_lst)]
                        if ms1_pr_df.shape[0] > 0:
                            # print('Number of MS1 pr mz in list:', ms1_pr_df.shape[0])
                            ms1_pr_df['ppm'] = abs(1e6 * (ms1_pr_df['mz'] - mz_lib) / mz_lib)
                            # select best intensity in the precursor ppm range. Priority: i > ppm
                            # ms1_pr_df = ms1_pr_df.sort_values(by=['i', 'ppm'], ascending=[False, True])
                            ms1_pr_df = ms1_pr_df.sort_values(by='i', ascending=False)
                            # print('ms1_pr_df')
                            # print(ms1_pr_df)
                            ms1_pr_se = ms1_pr_df.iloc[0]
                            ms1_mz = ms1_pr_se['mz']
                            ms1_i = ms1_pr_se['i']
                            ms1_pr_ppm = 1e6 * (ms1_mz - mz_lib) / mz_lib
                            # get spectra_df of corresponding MS2 DDA scan
                            if ms2_spec_idx in spectra_pl.items:
                                ms2_df = spectra_pl[ms2_spec_idx]
                                ms2_df = ms2_df.query('i > 0')
                                ms2_df = ms2_df.sort_values(by='i', ascending=False).reset_index(drop=True)
                                try:
                                    ms2_df.drop('rt', axis=1, inplace=True)
                                except (KeyError, ValueError):
                                    print('MS2_df do not have rt column...')
                            else:
                                print('!!!!!! MS2 spectra not in the list >>> >>>')
                        else:
                            print('!!!!!! Precursor m/z in MS1 not in the list >>> >>>')
                    else:
                        print('!!!!!! Precursor m/z in MS1 not in the list >>> >>>')
                else:
                    print('!!!!!! MS1 spectra not in the list >>> >>>')

                print('MS1 @ DDA#:%.0f | Total scan id:%.0f' % (ms2_dda_idx, ms1_spec_idx))
                print('MS2 @ DDA#:%.0f | Total scan id:%.0f' % (ms2_dda_idx, ms2_spec_idx))
                # print('--------------- NEXT _idx')

        # print('== == == == == == NEXT DF')

    else:
        print('=== ===DO NOT have this precursor pr_mz == %f and func_id == %f and scan_id == %f!!!!!!'
              % (mz, func_id, ms2_scan_id)
              )

    spec_info_dct = {'ms1_i': ms1_i, 'ms1_mz': ms1_mz, 'ms1_pr_ppm': ms1_pr_ppm, 'ms1_rt': ms1_rt, 'ms2_rt': ms2_rt,
                     '_ms1_spec_idx': ms1_spec_idx, '_ms2_spec_idx': ms2_spec_idx, 'ms1_df': ms1_df, 'ms2_df': ms2_df}

    return spec_info_dct


def get_xic_from_pl(xic_ms1_lst, ms1_xic_df, xic_ppm):

    ms1_xic_dct = {}

    xic_ms1_l_lst = ppm_window_para(xic_ms1_lst, -1 * xic_ppm)
    xic_ms1_h_lst = ppm_window_para(xic_ms1_lst, xic_ppm)
    xic_ms_info_lst = zip(xic_ms1_lst, xic_ms1_l_lst, xic_ms1_h_lst)

    for _xic_mz_info in xic_ms_info_lst:
        _xic_mz = _xic_mz_info[0]
        ms1_low = _xic_mz_info[1]
        ms1_high = _xic_mz_info[2]
        if _xic_mz > 0:
            ms1_query = '%f <= mz <= %f' % (ms1_low, ms1_high)
            # print(ms1_query)
            _found_ms1_df = ms1_xic_df.query(ms1_query)
            _found_ms1_df.loc[:, 'ppm'] = 1e6 * (_found_ms1_df['mz'] - _xic_mz) / _xic_mz
            _found_ms1_df.loc[:, 'ppm'] = _found_ms1_df['ppm'].abs()
            _found_ms1_df.loc[:, 'mz'] = _xic_mz
            _found_ms1_df.sort_values(by=['rt', 'i', 'ppm'],
                                      ascending=[True, False, True], inplace=True)
            _found_ms1_df.drop_duplicates(subset=['rt'], keep='first', inplace=True)
            # print('_found_ms1_df.shape', _found_ms1_df.shape)
            ms1_xic_dct[_xic_mz] = _found_ms1_df
    return ms1_xic_dct


def get_xic(ms1_mz, mzml, rt_range, ppm=500, ms1_precision=50e-6, msn_precision=500e-6, vendor='waters'):

    waters_obo_lst = [('MS:1000016', ['value']), ('MS:1000744', ['value']), ('MS:1000042', ['value']),
                      ('MS:1000796', ['value']), ('MS:1000514', ['name']), ('MS:1000515', ['name']),
                      ('MS:1000769', ['name']), ('MS:1000526', ['name']), ('MS:1002476', ['name']),
                      ('MS:1000514', ['name']), ('MS:1000515', ['name'])]
    thermo_obo_lst = [('MS:1000511', ['value']), ('MS:1000768', ['name']), ('MS:1000563', ['name'])]
    vendor_obo_lst = thermo_obo_lst + waters_obo_lst
    for _obo in vendor_obo_lst:
        if _obo not in pymzml.minimum.MIN_REQ:
            pymzml.minimum.MIN_REQ.append(_obo)
        else:
            pass
            # end hot patch

    rt_start = rt_range[0]
    rt_end = rt_range[1]

    ms1_low = ms1_mz - ms1_mz * ppm * 1e-6
    ms1_high = ms1_mz + ms1_mz * ppm * 1e-6
    ms1_query = '%f <= mz <= %f' % (ms1_low, ms1_high)
    print('Find XIC in range: %s' % ms1_query)

    ms1_xic_df = pd.DataFrame()

    spec_title_obo = 'MS:1000796'
    scan_rt_obo = 'MS:1000016'
    spec_level_obo = 'MS:1000511'

    spec_obj = pymzml.run.Reader(mzml, MS1_Precision=ms1_precision, MSn_Precision=msn_precision)

    if vendor == 'waters':

        scan_info_re = re.compile(r'(.*)(function=)(\d{1,2})(.*)(scan=)(\d*)(.*)')

        for _spectrum in spec_obj:

            if spec_title_obo in _spectrum.keys() and scan_rt_obo in _spectrum.keys():
                _spectrum_title = _spectrum[spec_title_obo]
                _scan_rt = float(_spectrum[scan_rt_obo])
                scan_info_checker = scan_info_re.match(_spectrum_title)

                if rt_start <= _scan_rt <= rt_end and scan_info_checker:
                    _function = int(scan_info_checker.groups()[2])
                    # _scan_id = int(scan_info_checker.groups()[5])
                    if _function == 1:
                        _tmp_spec_df = pd.DataFrame(data=_spectrum.peaks, columns=['mz', 'i'])
                        _found_ms1_df = _tmp_spec_df.query(ms1_query)
                        _found_ms1_df.loc[:, 'ppm'] = 1e6 * (_found_ms1_df['mz'] - ms1_mz) / ms1_mz
                        _found_ms1_df.loc[:, 'ppm'] = _found_ms1_df['ppm'].abs()
                        _found_ms1_df.loc[:, 'scan_time'] = _scan_rt

                        ms1_xic_df = ms1_xic_df.append(_found_ms1_df.sort_values(by='ppm').head(1))

    elif vendor == 'thermo':

        for _spectrum in spec_obj:

            if spec_level_obo in _spectrum.keys() and scan_rt_obo in _spectrum.keys():
                _spec_level = _spectrum[spec_level_obo]
                _scan_rt = float(_spectrum[scan_rt_obo])
                if rt_start <= _scan_rt <= rt_end:
                    if _spec_level == 1:
                        _tmp_spec_df = pd.DataFrame(data=_spectrum.peaks, columns=['mz', 'i'])
                        _found_ms1_df = _tmp_spec_df.query(ms1_query)
                        _found_ms1_df.loc[:, 'ppm'] = 1e6 * (_found_ms1_df['mz'] - ms1_mz) / ms1_mz
                        _found_ms1_df.loc[:, 'ppm'] = _found_ms1_df['ppm'].abs()
                        _found_ms1_df.loc[:, 'scan_time'] = _scan_rt

                        ms1_xic_df = ms1_xic_df.append(_found_ms1_df.sort_values(by='ppm').head(1))

    return ms1_xic_df


def get_xic_all(core_list, mzml, rt_range, ms1_precision=50e-6, msn_precision=500e-6, vendor='waters'):

    waters_obo_lst = [('MS:1000016', ['value']), ('MS:1000744', ['value']), ('MS:1000042', ['value']),
                      ('MS:1000796', ['value']), ('MS:1000514', ['name']), ('MS:1000515', ['name']),
                      ('MS:1000769', ['name']), ('MS:1000526', ['name'])]
    thermo_obo_lst = [('MS:1000511', ['value']), ('MS:1000768', ['name']), ('MS:1000563', ['name'])]
    vendor_obo_lst = thermo_obo_lst + waters_obo_lst
    for _obo in vendor_obo_lst:
        if _obo not in pymzml.minimum.MIN_REQ:
            pymzml.minimum.MIN_REQ.append(_obo)
        else:
            pass
            # end hot patch

    rt_start = rt_range[0]
    rt_end = rt_range[1]

    ms1_xic_dct = {}

    for _mz in core_list:

        ms1_xic_dct[_mz] = pd.DataFrame()
    spec_title_obo = 'MS:1000796'
    scan_rt_obo = 'MS:1000016'
    spec_level_obo = 'MS:1000511'

    ims_obo = 'MS:1002476'
    mz_array_obo = 'MS:1000514'
    i_array_obo = 'MS:1000515'

    spec_obj = pymzml.run.Reader(mzml, MS1_Precision=ms1_precision, MSn_Precision=msn_precision)

    if vendor == 'waters':
        try:
            scan_info_re = re.compile(r'(.*)(function=)(\d{1,2})(.*)')
            for _spectrum in spec_obj:
                if spec_title_obo in _spectrum.keys() and scan_rt_obo in _spectrum.keys():
                    _spectrum_title = _spectrum[spec_title_obo]
                    _scan_rt = float(_spectrum[scan_rt_obo])
                    scan_info_checker = scan_info_re.match(_spectrum_title)

                    if rt_start <= _scan_rt <= rt_end and scan_info_checker:
                        _function = int(scan_info_checker.groups()[2])
                        if _function == 1:
                            print('Reading MS survey scan @:', _scan_rt)
                            # slow but more accurate mode. At least 10 time slower
                            # _tmp_spec_df = pd.DataFrame(data=_spectrum.peaks, columns=['mz', 'i'])
                            for _ms1_xic in core_list:

                                ms1_xic_df = ms1_xic_dct[_ms1_xic]
                                _tmp_ms1_xic_df = ms1_xic_df.copy()

                                # faster mode
                                _xic_lst = _spectrum.hasPeak(_ms1_xic)
                                if len(_xic_lst) == 1:
                                    _tmp_mz_df = pd.DataFrame(data=_xic_lst, columns=['mz', 'i'])
                                    _tmp_mz_df.loc[:, 'rt'] = _scan_rt
                                    _tmp_mz_df.loc[:, 'mz'] = _ms1_xic
                                    ms1_xic_df = _tmp_ms1_xic_df.append(_tmp_mz_df)
                                    ms1_xic_dct[_ms1_xic] = ms1_xic_df

                                if len(_xic_lst) > 1:
                                    _tmp_mz_df = pd.DataFrame(data=_xic_lst, columns=['mz', 'i'])
                                    _tmp_mz_df.loc[:, 'rt'] = _scan_rt
                                    _tmp_mz_df.loc[:, 'mz'] = _ms1_xic
                                    ms1_xic_df = _tmp_ms1_xic_df.append(_tmp_mz_df.sort_values(by='i',
                                                                                               ascending=False).head(1))
                                    ms1_xic_dct[_ms1_xic] = ms1_xic_df
        except KeyError:
            print('!! Spectra format does not match to vendor settings !! Please check your file and settings !!')

    elif vendor == 'thermo':
        print('Thermo files')
        try:
            for _spectrum in spec_obj:

                if spec_level_obo in _spectrum.keys() and scan_rt_obo in _spectrum.keys():
                    # ms_level = _spectrum[spec_level_obo]
                    _spectrum_level = _spectrum[spec_level_obo]
                    _scan_rt = float(_spectrum[scan_rt_obo])

                    if rt_start <= _scan_rt <= rt_end:
                        if _spectrum_level == 1:
                            print('Reading MS survey scan @:', _scan_rt)
                            # slow but more accurate mode. At least 10 time slower
                            # _tmp_spec_df = pd.DataFrame(data=_spectrum.peaks, columns=['mz', 'i'])
                            for _ms1_xic in core_list:

                                ms1_xic_df = ms1_xic_dct[_ms1_xic]
                                _tmp_ms1_xic_df = ms1_xic_df.copy()

                                # faster mode
                                _xic_lst = _spectrum.hasPeak(_ms1_xic)
                                if len(_xic_lst) == 1:
                                    _tmp_mz_df = pd.DataFrame(data=_xic_lst, columns=['mz', 'i'])
                                    _tmp_mz_df.loc[:, 'rt'] = _scan_rt
                                    _tmp_mz_df.loc[:, 'mz'] = _ms1_xic
                                    ms1_xic_df = _tmp_ms1_xic_df.append(_tmp_mz_df)
                                    ms1_xic_dct[_ms1_xic] = ms1_xic_df

                                if len(_xic_lst) > 1:
                                    _tmp_mz_df = pd.DataFrame(data=_xic_lst, columns=['mz', 'i'])
                                    _tmp_mz_df.loc[:, 'rt'] = _scan_rt
                                    _tmp_mz_df.loc[:, 'mz'] = _ms1_xic
                                    ms1_xic_df = _tmp_ms1_xic_df.append(_tmp_mz_df.sort_values
                                                                        (by='i', ascending=False).head(1))
                                    ms1_xic_dct[_ms1_xic] = ms1_xic_df
        except KeyError:
            print('!! Spectra format does not match to vendor settings !! Please check your file and settings !!')

    return ms1_xic_dct


def get_spec_info(lpp_all_group_key_lst, checked_info_groups, usr_scan_info_df):
    lpp_spec_info_dct = {}
    for group_key in lpp_all_group_key_lst:
        _subgroup_df = checked_info_groups.get_group(group_key)
        _samemz_se = _subgroup_df.iloc[0, :].squeeze()
        _usr_ms2_pr_mz = _samemz_se['MS2_PR_mz']

        _usr_ms2_dda_rank = _samemz_se['DDA_rank']
        _usr_ms2_scan_id = _samemz_se['scan_number']
        _usr_mz_lib = _samemz_se['Lib_mz']
        _tmp_chk_df = usr_scan_info_df.query('MS2_PR_mz == %.6f and DDA_rank == %i and scan_number == %i'
                                             % (_usr_ms2_pr_mz, _usr_ms2_dda_rank, _usr_ms2_scan_id))
        if _tmp_chk_df.shape[0] == 1:
            _tmp_info_dct = {'MS2_PR_mz': _usr_ms2_pr_mz, 'DDA_rank': _usr_ms2_dda_rank,
                             'scan_number': _usr_ms2_scan_id, 'Lib_mz': _usr_mz_lib}
            lpp_spec_info_dct[group_key] = _tmp_info_dct

    return lpp_spec_info_dct

# if __name__ == '__main__':
#
#     usr_mzml = r'test\CM_neg_30min.mzML'
#     usr_dda_top = 12
#     usr_rt_range = [25, 27]
#
#     usr_scan_info_df, usr_spec_pl = extract_mzml(usr_mzml, usr_rt_range, usr_dda_top)
#
#     print(usr_scan_info_df.head(5))
#     print(usr_spec_pl.items)
