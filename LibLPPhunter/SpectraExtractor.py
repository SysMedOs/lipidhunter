# -*- coding: utf-8 -*-
# Copyright 2015-2016 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

from __future__ import division
from __future__ import print_function
import re

import pandas as pd
import pymzml


def extract_mzml(mzml, rt_range, dda_top=6, ms1_precision=50e-6, msn_precision=500e-6):

    """
    Extract mzML to a scan info DataFrame and a pandas panel for spectra DataFrame of mz and i

    :param mzml: The file path of mzML file
    :type mzml: str
    :param rt_range: A List of RT. e.g. [15, 30] for 15 to 30 min
    :type rt_range: list
    :param dda_top: DDA settings e.g. DDA TOP 6
    :type dda_top: int
    :param ms1_precision: e.g. 50e-6 for 50 ppm
    :type ms1_precision: float
    :param msn_precision: e.g. 500e-6 for 500 ppm
    :type msn_precision: float
    :returns: scan_info_df, spec_pl
    :rtype: pandas.DataFrame, pandas.Panel

    """

    waters_obo_lst = (('MS:1000016', ['value']), ('MS:1000744', ['value']), ('MS:1000042', ['value']),
                      ('MS:1000796', ['value']), ('MS:1000514', ['name']), ('MS:1000515', ['name']))
    for _obo in waters_obo_lst:
        if _obo not in pymzml.minimum.MIN_REQ:
            pymzml.minimum.MIN_REQ.append(_obo)
        else:
            pass

    rt_start = rt_range[0]
    rt_end = rt_range[1]

    print('==> Start to process file: %s' % mzml)
    print('=== ==> RT: %.2f -> %.2f with DDA Top % i' % (rt_start, rt_end, dda_top))

    spec_title_obo = 'MS:1000796'
    scan_rt_obo = 'MS:1000016'
    scan_pr_mz_obo = 'MS:1000744'

    spec_obj = pymzml.run.Reader(mzml, MS1_Precision=ms1_precision, MSn_Precision=msn_precision)

    scan_info_re = re.compile(r'(.*)(function=)(\d{1,2})(.*)(scan=)(\d*)(.*)')

    spec_idx = 0
    dda_event_idx = 0
    spec_idx_lst = []
    dda_event_lst = []
    rt_lst = []
    function_lst = []
    scan_id_lst = []
    pr_mz_lst = []

    scan_info_dct = {'spec_index': spec_idx_lst, 'rt': rt_lst, 'dda_event_idx': dda_event_lst,
                     'function': function_lst, 'scan_id': scan_id_lst, 'pr_mz': pr_mz_lst}

    spec_dct = {}

    ms2_function_range_lst = range(2, dda_top + 1)
    function_range_lst = range(1, dda_top + 1)

    for _spectrum in spec_obj:
        pr_mz = 0
        if spec_title_obo in _spectrum.keys() and scan_rt_obo in _spectrum.keys():
            _spectrum_title = _spectrum[spec_title_obo]
            _scan_rt = float(_spectrum[scan_rt_obo])
            scan_info_checker = scan_info_re.match(_spectrum_title)

            if rt_start <= _scan_rt <= rt_end and scan_info_checker:
                _function = int(scan_info_checker.groups()[2])
                _scanid = int(scan_info_checker.groups()[5])
                if _function in function_range_lst:
                    _tmp_spec_df = pd.DataFrame(data=_spectrum.peaks, columns=['mz', 'i'])
                    if _function == 1:
                        dda_event_idx += 1

                        _tmp_spec_df = _tmp_spec_df.sort_values(by='i', ascending=False).head(1000)
                        _tmp_spec_df = _tmp_spec_df.reset_index(drop=True)

                        spec_dct[spec_idx] = _tmp_spec_df

                    if _function in ms2_function_range_lst:
                        pr_mz = _spectrum[scan_pr_mz_obo]
                        spec_dct[spec_idx] = _tmp_spec_df.query('i >= 10')

                    spec_idx_lst.append(spec_idx)
                    dda_event_lst.append(dda_event_idx)
                    rt_lst.append(_scan_rt)
                    function_lst.append(_function)
                    scan_id_lst.append(_scanid)
                    pr_mz_lst.append(pr_mz)

                    # print(dda_event_idx, spec_idx, _ms2_function, _scanid, _scan_rt, pr_mz)

                    spec_idx += 1

    scan_info_df = pd.DataFrame(data=scan_info_dct, columns=['dda_event_idx', 'spec_index', 'rt',
                                                             'function', 'scan_id', 'pr_mz'])

    scan_info_df = scan_info_df.sort_values(by='rt')
    scan_info_df = scan_info_df.round({'pr_mz': 6})

    # scan_info_df.to_excel('scan_info_df.xlsx')

    spec_pl = pd.Panel(data=spec_dct)
    print('=== ==> --> mzML extracted')

    return scan_info_df, spec_pl


def get_spectra(mz, function, ms2_scan_id, scan_info_df, spectra_pl, dda_top=12):

    ms1_df = pd.DataFrame()
    ms2_df = pd.DataFrame()
    ms1_spec_idx = 0
    ms2_spec_idx = 0
    ms1_rt = 0
    ms2_rt = 0

    if mz in scan_info_df['pr_mz'].tolist():
        _tmp_mz_scan_info_df = scan_info_df.query('pr_mz == %f and function == %f and scan_id == %f'
                                                  % (mz, function, ms2_scan_id)
                                                  )

        if _tmp_mz_scan_info_df.shape[0] == 1:
            ms2_spec_idx = _tmp_mz_scan_info_df.get_value(_tmp_mz_scan_info_df.index[0], 'spec_index')
            ms2_dda_idx = _tmp_mz_scan_info_df.get_value(_tmp_mz_scan_info_df.index[0], 'dda_event_idx')
            ms2_function = _tmp_mz_scan_info_df.get_value(_tmp_mz_scan_info_df.index[0], 'function')
            ms2_scan_id = _tmp_mz_scan_info_df.get_value(_tmp_mz_scan_info_df.index[0], 'scan_id')
            ms2_rt = _tmp_mz_scan_info_df.get_value(_tmp_mz_scan_info_df.index[0], 'rt')

            print('%.6f @ DDA#:%.0f | Total scan id:%.0f | function: %.0f | Scan ID: %.0f | RT: %.4f'
                  % (mz, ms2_dda_idx, ms2_spec_idx, ms2_function, ms2_scan_id, ms2_rt)
                  )

            # get spectra_df of corresponding MS survey scan
            tmp_ms1_info_df = scan_info_df.query('dda_event_idx == %i and function == 1' % ms2_dda_idx)
            ms1_spec_idx = tmp_ms1_info_df['spec_index'].tolist()[0]
            ms1_rt = tmp_ms1_info_df['rt'].tolist()[0]
            if ms1_spec_idx in spectra_pl.items:
                ms1_df = spectra_pl[ms1_spec_idx]
                ms1_df = ms1_df.query('i > 0')
                ms1_df = ms1_df.sort_values(by='i', ascending=False).reset_index(drop=True)
            else:
                pass

            # get spectra_df of corresponding MS2 DDA scan
            if ms2_spec_idx in spectra_pl.items:
                ms2_df = spectra_pl[ms2_spec_idx]
                ms2_df = ms2_df.query('i > 0')
                ms2_df = ms2_df.sort_values(by='i', ascending=False).reset_index(drop=True)
            else:
                pass

            # if _ms1_spec_idx - _ms2_spec_idx > dda_top:
            #     print('!!!!!!!!!!!! MS1 is NOT for this MS/MS !!!!!!!!!!!!')

            print('MS1 @ DDA#:%.0f | Total scan id:%.0f' % (ms2_dda_idx, ms1_spec_idx))
            # print(ms1_df.head(5))
            print('MS2 @ DDA#:%.0f | Total scan id:%.0f' % (ms2_dda_idx, ms2_spec_idx))
            # print(ms2_df.head(5))

            print('--------------- NEXT _idx')
        print('== == == == == == NEXT DF')

    else:
        print('=== ===DO NOT have this precursor pr_mz == %f and function == %f and scan_id == %f!!!!!!'
              % (mz, function, ms2_scan_id)
              )
    print(ms2_rt)
    return ms1_rt, ms2_rt, ms1_spec_idx, ms2_spec_idx, ms1_df, ms2_df


def get_xic(ms1_mz, mzml, rt_range, ppm=500, ms1_precision=50e-6, msn_precision=500e-6):

    waters_obo_lst = (('MS:1000016', ['value']), ('MS:1000744', ['value']), ('MS:1000042', ['value']),
                      ('MS:1000796', ['value']), ('MS:1000514', ['name']), ('MS:1000515', ['name']))
    for _obo in waters_obo_lst:
        if _obo not in pymzml.minimum.MIN_REQ:
            pymzml.minimum.MIN_REQ.append(_obo)
        else:
            pass

    rt_start = rt_range[0]
    rt_end = rt_range[1]

    ms1_low = ms1_mz - ms1_mz * ppm * 1e-6
    ms1_high = ms1_mz + ms1_mz * ppm * 1e-6
    ms1_query = '%f <= mz <= %f' % (ms1_low, ms1_high)
    print('Find XIC in range: %s' % ms1_query)

    ms1_xic_df = pd.DataFrame()

    spec_title_obo = 'MS:1000796'
    scan_rt_obo = 'MS:1000016'

    spec_obj = pymzml.run.Reader(mzml, MS1_Precision=ms1_precision, MSn_Precision=msn_precision)

    scan_info_re = re.compile(r'(.*)(function=)(\d{1,2})(.*)(scan=)(\d*)(.*)')

    for _spectrum in spec_obj:

        if spec_title_obo in _spectrum.keys() and scan_rt_obo in _spectrum.keys():
            _spectrum_title = _spectrum[spec_title_obo]
            _scan_rt = float(_spectrum[scan_rt_obo])
            scan_info_checker = scan_info_re.match(_spectrum_title)

            if rt_start <= _scan_rt <= rt_end and scan_info_checker:
                _function = int(scan_info_checker.groups()[2])
                _scanid = int(scan_info_checker.groups()[5])
                if _function == 1:
                    _tmp_spec_df = pd.DataFrame(data=_spectrum.peaks, columns=['mz', 'i'])
                    _found_ms1_df = _tmp_spec_df.query(ms1_query)
                    _found_ms1_df['ppm'] = 1e6 * (_found_ms1_df['mz'] - ms1_mz) / ms1_mz
                    _found_ms1_df['ppm'] = _found_ms1_df['ppm'].abs()
                    _found_ms1_df['rt'] = _scan_rt

                    ms1_xic_df = ms1_xic_df.append(_found_ms1_df.sort_values(by='ppm').head(1))

    return ms1_xic_df


if __name__ == '__main__':

    usr_mzml = r'D:\project_mzML\CM_DDA_neg_mzML\070120_CM_neg_30min_SIN_II.mzML'
    usr_dda_top = 12
    usr_rt_range = [25, 27]

    usr_scan_info_df, usr_spec_pl = extract_mzml(usr_mzml, usr_rt_range, usr_dda_top)

    print(usr_scan_info_df.head(5))
    print(usr_spec_pl.items)
