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

import getopt
import math
import os
import sys
import time
from multiprocessing import Pool

import pandas as pd

from LibLipidHunter.LipidComposer import LipidComposer
from LibLipidHunter.SpectraReader import extract_mzml
from LibLipidHunter.SpectraReader import get_spectra
from LibLipidHunter.SpectraReader import get_xic_from_pl
from LibLipidHunter.SpectraReader import get_spec_info
from LibLipidHunter.LogPageCreator import LogPageCreator
from LibLipidHunter.PrecursorHunter import PrecursorHunter
from LibLipidHunter.ScoreHunter import get_lipid_info


def huntlipids(param_dct, error_lst):
    """

    :param param_dct:
    example
    hunter_param_dct = {'fawhitelist_path_str': r'D:\lipidhunter\ConfigurationFiles\FA_Whitelist.xlsx',
                        'mzml_path_str': r'D:\lipidhunter\mzML\MS2\070120_CM_neg_70min_SIN_I.mzML',
                        'img_output_folder_str': r'D:\lipidhunter\Temp\Test2',
                        'xlsx_output_path_str': r'D:\lipidhunter\Temp\Test2\t2.xlsx',
                        'lipid_specific_cfg': r'D:\lipidhunter\ConfigurationFiles\PL_specific_ion_cfg.xlsx',
                        'hunter_start_time': '2017-12-21_15-27-49',
                        'vendor': 'waters', 'experiment_mode': 'LC-MS', 'lipid_type': 'PC', 'charge_mode': '[M+HCOO]-',
                        'rt_start': 20.0, 'rt_end': 25.0, 'mz_start': 700.0, 'mz_end': 800.0,
                        'rank_score': True, 'rank_score_filter': 27.5, 'score_filter': 27.5,
                        'isotope_score_filter': 75.0, 'fast_isotope': False,
                        'ms_th': 1000, 'ms_ppm': 20, 'ms_max': 0, 'pr_window': 0.75, 'dda_top': 6,
                        'ms2_th': 10, 'ms2_ppm': 50, 'ms2_infopeak_threshold': 0.001,
                        'hg_th': 10.0, 'hg_ppm': 100.0, 'ms2_hginfopeak_threshold': 0.001,
                        'score_cfg': r'D:\lipidhunter\ConfigurationFiles\Score_cfg.xlsx',
                        'fa_white_list_cfg': r'D:\lipidhunter\ConfigurationFiles\FA_Whitelist.xlsx',
                        'hunter_folder': r'D:\lipidhunter',
                        'core_number': 3, 'max_ram': 5, 'img_type': u'png', 'img_dpi': 300}
    :param error_lst: empty list to store error info to display on GUI.
    :return:
    """

    start_time = time.clock()
    lipidcomposer = LipidComposer()

    usr_lipid_class = param_dct['lipid_type']
    usr_charge = param_dct['charge_mode']
    usr_vendor = param_dct['vendor']
    usr_fa_xlsx = param_dct['fawhitelist_path_str']
    usr_mzml = param_dct['mzml_path_str']
    output_folder = param_dct['img_output_folder_str']
    output_sum_xlsx = param_dct['xlsx_output_path_str']

    key_frag_cfg = param_dct['lipid_specific_cfg']
    score_cfg = param_dct['score_cfg']

    usr_rt_range = [param_dct['rt_start'], param_dct['rt_end']]
    # usr_pr_mz_range = [param_dct['mz_start'], param_dct['mz_end']]
    mz_start = param_dct['mz_start']
    mz_end = param_dct['mz_end']
    usr_dda_top = param_dct['dda_top']
    usr_ms1_threshold = param_dct['ms_th']
    usr_ms1_max = param_dct['ms_max']
    usr_ms2_threshold = param_dct['ms2_th']
    # usr_ms2_hg_threshold = param_dct['hg_th']
    usr_ms1_ppm = param_dct['ms_ppm']
    usr_ms2_ppm = param_dct['ms2_ppm']
    usr_ms1_precision = usr_ms1_ppm * 1e-6
    usr_ms2_precision = usr_ms2_ppm * 1e-6
    # usr_ms2_hg_precision = param_dct['hg_ppm'] * 1e-6
    # usr_rank_score_filter = param_dct['rank_score_filter']
    # usr_score_filter = param_dct['score_filter']
    # usr_isotope_score_filter = param_dct['isotope_score_filter']
    # usr_ms2_info_th = param_dct['ms2_infopeak_threshold']
    # usr_ms2_hginfo_th = param_dct['ms2_hginfopeak_threshold']
    # usr_rank_mode = param_dct['rank_score']
    # usr_fast_isotope = param_dct['fast_isotope']

    hunter_start_time_str = param_dct['hunter_start_time']

    # keep stay in current working directory
    current_path = os.getcwd()
    if os.path.isdir(output_folder):
        os.chdir(output_folder)
        if os.path.isdir('LipidHunter_Results_Figures_%s' % hunter_start_time_str):
            pass
        else:
            os.mkdir('LipidHunter_Results_Figures_%s' % hunter_start_time_str)
    os.chdir(current_path)

    print('=== ==> --> Start to process >>>')
    print('=== ==> --> Phospholipid class: %s >>>' % usr_lipid_class)

    composer_param_dct = {'fa_whitelist': usr_fa_xlsx, 'lipid_type': usr_lipid_class,
                          'charge_mode': usr_charge, 'exact_position': 'FALSE'}
    usr_lipid_master_df = lipidcomposer.compose_lipid(param_dct=composer_param_dct)
    usr_fa_df = lipidcomposer.calc_fa_query(usr_lipid_class, usr_fa_xlsx, ms2_ppm=usr_ms2_ppm)

    print('=== ==> --> Lipid Master table generated >>>', usr_lipid_master_df.shape)
    # print(usr_lipid_master_df.head())

    # parameters from settings tab
    usr_core_num = param_dct['core_number']
    usr_max_ram = param_dct['max_ram']

    lipid_info_df = usr_lipid_master_df

    # cut lib info to the user defined m/z range
    pos_charge_lst = ['[M+H]+', '[M+Na]+', '[M+NH4]+']
    neg_charge_lst = ['[M-H]-', '[M+HCOO]-', '[M+CH3COO]-']
    if usr_charge in neg_charge_lst:
        if usr_lipid_class == 'PC':
            if usr_charge == '[M+HCOO]-':
                lipid_info_df = lipid_info_df[(mz_start <= lipid_info_df['[M+HCOO]-_MZ'])
                                              & (lipid_info_df['[M+HCOO]-_MZ'] <= mz_end)]
            elif usr_charge == '[M+CH3COO]-':
                lipid_info_df = lipid_info_df[(mz_start <= lipid_info_df['[M+CH3COO]-_MZ'])
                                              & (lipid_info_df['[M+CH3COO]-_MZ'] <= mz_end)]
            else:
                error_lst.append('PC charge not supported.  User input charge = %s. '
                                 'LipidHunter support [M+HCOO]- and [M+CH3COO]-.' % usr_charge)
        else:
            lipid_info_df = lipid_info_df[
                (mz_start <= lipid_info_df['[M-H]-_MZ']) & (lipid_info_df['[M-H]-_MZ'] <= mz_end)]
    elif usr_charge in pos_charge_lst:
        if usr_lipid_class == 'TG':
            # TODO(zhixu.ni@uni-leipzig.de): @Georgia add TG here please :)
            pass
    else:
        error_lst.append('Lipid class or charge NOT supported.  User input lipid class = %s, charge = %s. '
                         % (usr_lipid_class, usr_charge))

    # TODO(zhixu.ni@uni-leipzig.de): Add more error to the error_lst.

    pr_hunter = PrecursorHunter(lipid_info_df, param_dct)

    # keep stay in current working directory
    current_path = os.getcwd()
    if os.path.isdir(output_folder):
        os.chdir(output_folder)
        if os.path.isdir('LipidHunter_Results_Figures_%s' % hunter_start_time_str):
            print('... Output folder existed...')
        else:
            os.mkdir('LipidHunter_Results_Figures_%s' % hunter_start_time_str)
            print('... Output folder created...')
    else:
        os.mkdir(output_folder)
        os.chdir(output_folder)
        os.mkdir('LipidHunter_Results_Figures_%s' % hunter_start_time_str)
        print('... Output folder created...')
    os.chdir(current_path)

    log_pager = LogPageCreator(output_folder, hunter_start_time_str, param_dct)

    output_df = pd.DataFrame()

    print('=== ==> --> Start to process')
    print('=== ==> --> Phospholipid class: %s' % usr_lipid_class)

    # generate the indicator table

    usr_weight_df = pd.read_excel(score_cfg, index_col='Type')

    usr_key_frag_df = pd.read_excel(key_frag_cfg)
    usr_key_frag_df = usr_key_frag_df.query('EXACTMASS > 0')

    # get the information from the following columns and leave the rest
    usr_key_frag_df = usr_key_frag_df[['CLASS', 'TYPE', 'EXACTMASS', 'PR_CHARGE', 'LABEL', 'CHARGE_MODE']]

    print('=== ==> --> Start to parse mzML')
    # extract all spectra from mzML to pandas DataFrame
    usr_scan_info_df, usr_spectra_pl, ms1_xic_df = extract_mzml(usr_mzml, usr_rt_range, dda_top=usr_dda_top,
                                                                ms1_threshold=usr_ms1_threshold,
                                                                ms2_threshold=usr_ms2_threshold,
                                                                ms1_precision=usr_ms1_precision,
                                                                ms2_precision=usr_ms2_precision,
                                                                vendor=usr_vendor, ms1_max=usr_ms1_max)

    print('MS1_XIC_df.shape', ms1_xic_df.shape)

    ms1_obs_pr_df, sub_pl_group_lst = pr_hunter.get_matched_pr(usr_scan_info_df, usr_spectra_pl, ms1_max=usr_ms1_max,
                                                               core_num=usr_core_num, max_ram=usr_max_ram)

    if isinstance(ms1_obs_pr_df, str):
        return '!! NO suitable precursor --> Check settings!!\n'

    print('=== ==> --> ms1 precursor matched')

    # remove bad precursors
    checked_info_df = pd.DataFrame()
    for _idx, _check_scan_se in usr_scan_info_df.iterrows():
        _dda_rank = _check_scan_se['DDA_rank']
        _scan_id = _check_scan_se['scan_number']
        _tmp_usr_df = ms1_obs_pr_df.query('DDA_rank == %f and scan_number == %f' % (_dda_rank, _scan_id))
        checked_info_df = checked_info_df.append(_tmp_usr_df)

    checked_info_df.sort_values(by=['MS2_PR_mz', 'scan_number'])

    ms1_xic_mz_lst = ms1_obs_pr_df['MS1_XIC_mz'].tolist()
    ms1_xic_mz_lst = sorted(set(ms1_xic_mz_lst))
    print('ms1_xic_mz_lst', len(ms1_xic_mz_lst))
    print(ms1_xic_mz_lst)

    print('=== ==> --> Start to extract XIC')
    if len(ms1_xic_mz_lst) >= 3 * usr_core_num:
        sub_len = int(math.ceil(len(ms1_xic_mz_lst) / usr_core_num))
        core_key_list = map(None, *(iter(ms1_xic_mz_lst),) * sub_len)
    else:
        core_key_list = [ms1_xic_mz_lst]

    # Start multiprocessing to get XIC
    print('!!!!!! Start multiprocessing to get XIC ==> ==> ==> Number of Cores: %i' % usr_core_num)
    xic_dct = {}

    if 1 < usr_core_num < len(core_key_list):
        parallel_pool = Pool(usr_core_num)
        xic_results_lst = []
        core_worker_count = 1
        for core_list in core_key_list:
            if isinstance(core_list, tuple) or isinstance(core_list, list):
                if None in core_list:
                    core_list = filter(lambda x: x is not None, core_list)
                else:
                    pass
                print('>>> >>> Core #%i ==> ...... processing ......' % core_worker_count)
                print(core_list)
                xic_result = parallel_pool.apply_async(get_xic_from_pl, args=(core_list, ms1_xic_df, 500))
                core_worker_count += 1
                xic_results_lst.append(xic_result)

        parallel_pool.close()
        parallel_pool.join()

        for xic_result in xic_results_lst:
            try:
                sub_xic_dct = xic_result.get()
                if len(sub_xic_dct.keys()) > 0:
                    xic_dct = dict(xic_dct, **sub_xic_dct)
            except (KeyError, SystemError, ValueError):
                pass
    else:
        print('Using single core mode...')
        core_worker_count = 1
        for core_list in core_key_list:
            if isinstance(core_list, tuple) or isinstance(core_list, list):
                if None in core_list:
                    core_list = filter(lambda x: x is not None, core_list)
                else:
                    pass
                print('>>> >>> Core #%i ==> ...... processing ......' % core_worker_count)
                print(core_list)
                sub_xic_dct = get_xic_from_pl(core_list, ms1_xic_df, 500)
                core_worker_count += 1
                if len(sub_xic_dct.keys()) > 0:
                    xic_dct = dict(xic_dct, **sub_xic_dct)

    print('xic_dct', len(xic_dct.keys()))
    print(xic_dct.keys())

    if len(xic_dct.keys()) == 0:
        print('No precursor for XIC found')
        return '!! NO suitable precursor --> Check settings!!\n'
    else:
        print('=== ==> --> Number of XIC extracted: %i' % len(xic_dct.keys()))

    target_ident_lst = []
    checked_info_df.sort_values(by=['Lib_mz', 'scan_time', 'MS2_PR_mz'], ascending=[True, True, True], inplace=True)

    print('=== ==> --> Start to Hunt for Lipids !!')
    checked_info_groups = checked_info_df.groupby(['Lib_mz', 'MS2_PR_mz', 'Formula', 'scan_time', 'Ion'])
    lipid_all_group_key_lst = checked_info_groups.groups.keys()
    # lipid_all_group_key_lst = sorted(lipid_all_group_key_lst, key=lambda x: x[0])

    spec_sub_len = int(math.ceil(len(lipid_all_group_key_lst) / usr_core_num))
    spec_sub_key_lst = map(None, *(iter(lipid_all_group_key_lst),) * spec_sub_len)

    lipid_spec_info_dct = {}
    # usr_core_num = 1
    if usr_core_num > 1:
        parallel_pool = Pool(usr_core_num)
        spec_results_lst = []
        core_worker_count = 1
        for _sub_lst in spec_sub_key_lst:
            if isinstance(_sub_lst, tuple) or isinstance(_sub_lst, list):
                if None in _sub_lst:
                    _sub_lst = filter(lambda x: x is not None, _sub_lst)
                else:
                    pass
                print('>>> >>> Core #%i ==> ...... processing ......' % core_worker_count)
                spec_result = parallel_pool.apply_async(get_spec_info, args=(_sub_lst, checked_info_groups,
                                                                             usr_scan_info_df))
                core_worker_count += 1
                spec_results_lst.append(spec_result)

        parallel_pool.close()
        parallel_pool.join()

        for spec_result in spec_results_lst:
            try:
                sub_spec_dct = spec_result.get()
                if len(sub_spec_dct.keys()) > 0:
                    lipid_spec_info_dct = dict(lipid_spec_info_dct, **sub_spec_dct)
            except (KeyError, SystemError, ValueError):
                print('ValueError: must supply a tuple to get_group with multiple grouping keys')
    else:
        print('Using single core mode...')
        core_worker_count = 1
        for _sub_lst in spec_sub_key_lst:
            if isinstance(_sub_lst, tuple) or isinstance(_sub_lst, list):
                if None in _sub_lst:
                    _sub_lst = filter(lambda x: x is not None, _sub_lst)
                else:
                    pass
                print('>>> >>> Core #%i ==> ...... processing ......' % core_worker_count)
                sub_spec_dct = get_spec_info(_sub_lst, checked_info_groups, usr_scan_info_df)
                core_worker_count += 1
                if len(sub_spec_dct.keys()) > 0:
                    lipid_spec_info_dct = dict(lipid_spec_info_dct, **sub_spec_dct)

    print('lipid_spec_info_dct', len(lipid_spec_info_dct.keys()))

    # Single process ONLY. usr_spectra_pl is too big in RAM --> RAM leaking during copy
    lipid_spec_dct = {}
    spec_info_key_lst = lipid_spec_info_dct.keys()
    for _spec_group_key in spec_info_key_lst:
        _spec_info_dct = lipid_spec_info_dct[_spec_group_key]
        _usr_ms2_pr_mz = _spec_info_dct['MS2_PR_mz']
        _usr_ms2_dda_rank = _spec_info_dct['DDA_rank']
        _usr_ms2_scan_id = _spec_info_dct['scan_number']
        _usr_mz_lib = _spec_info_dct['Lib_mz']
        # TODO(zhixu.ni@uni-leipzig.de): change get_spectra into multiple processing mode
        usr_spec_info_dct = get_spectra(_usr_ms2_pr_mz, _usr_mz_lib, _usr_ms2_dda_rank, _usr_ms2_scan_id,
                                        ms1_xic_mz_lst, usr_scan_info_df, usr_spectra_pl,
                                        dda_top=usr_dda_top, ms1_precision=usr_ms1_precision, vendor=usr_vendor)
        lipid_spec_dct[_spec_group_key] = usr_spec_info_dct

    found_spec_key_lst = lipid_spec_dct.keys()
    found_spec_key_lst = sorted(found_spec_key_lst, key=lambda x: x[0])
    spec_key_num = len(found_spec_key_lst)
    lipid_part_key_lst = []
    if spec_key_num > (usr_core_num * 40):
        lipid_part_len = int(math.ceil(spec_key_num / 8))
        lipid_part_lst = map(None, *(iter(found_spec_key_lst),) * lipid_part_len)
        for part_lst in lipid_part_lst:
            if None in part_lst:
                part_lst = filter(lambda x: x is not None, part_lst)
            lipid_sub_len = int(math.ceil(len(part_lst) / usr_core_num))
            lipid_sub_key_lst = map(None, *(iter(part_lst),) * lipid_sub_len)
            lipid_part_key_lst.append(lipid_sub_key_lst)

    else:
        lipid_sub_len = int(math.ceil(spec_key_num / usr_core_num))
        lipid_sub_key_lst = map(None, *(iter(found_spec_key_lst),) * lipid_sub_len)
        lipid_part_key_lst.append(lipid_sub_key_lst)

    part_tot = len(lipid_part_key_lst)
    part_counter = 1

    # parse specific peak info
    pl_class_lst = ['PA', 'PC', 'PE', 'PG', 'PI', 'PS', 'PIP']
    pl_neg_chg_lst = ['[M-H]-', '[M+HCOO]-', '[M+CH3COO]-']
    if usr_lipid_class in pl_class_lst and usr_charge in pl_neg_chg_lst:
        charge_mode = 'NEG'
        usr_key_frag_df = pd.read_excel(key_frag_cfg)
        usr_key_frag_df = usr_key_frag_df.query('EXACTMASS > 0')
        # get the information from the following columns
        usr_key_frag_df = usr_key_frag_df[['CLASS', 'TYPE', 'EXACTMASS', 'PR_CHARGE', 'LABEL', 'CHARGE_MODE']]
        # find key peaks for the target PL class
        target_frag_df = usr_key_frag_df.query(r'CLASS == "%s" and TYPE == "FRAG" and PR_CHARGE == "%s"'
                                               % (usr_lipid_class, usr_charge))
        target_nl_df = usr_key_frag_df.query(r'CLASS == "%s" and TYPE == "NL" and PR_CHARGE == "%s"'
                                             % (usr_lipid_class, usr_charge))
        # add precursor to the list
        target_pr_df = pd.DataFrame(data={'CLASS': usr_lipid_class, 'TYPE': 'NL', 'EXACTMASS': 0.0,
                                          'PR_CHARGE': usr_charge, 'LABEL': 'PR', 'CHARGE_MODE': 'NEG'}, index=['PR'])
        target_nl_df = target_nl_df.append(target_pr_df)
        target_nl_df.reset_index(drop=True, inplace=True)

        # extract info for other classes
        other_frag_df = usr_key_frag_df.query('CLASS != "%s" and TYPE == "FRAG" and CHARGE_MODE == "%s"'
                                              % (usr_lipid_class, charge_mode))
        other_nl_df = usr_key_frag_df.query('CLASS != "%s" and TYPE == "NL" and CHARGE_MODE == "%s"'
                                            % (usr_lipid_class, charge_mode))
        key_frag_dct = {'target_frag_df': target_frag_df, 'target_nl_df': target_nl_df,
                        'other_frag_df': other_frag_df, 'other_nl_df': other_nl_df}
    else:
        key_frag_dct = {}

    print('... Key FRAG Dict Generated ...')

    for lipid_sub_key_lst in lipid_part_key_lst:

        if part_tot == 1:
            print('>>> Start multiprocessing to get Score ==> Max Number of Cores: %i' % usr_core_num)
        else:
            print('>>> Start multiprocessing to get Score ==> Part %i / %i --> Max Number of Cores: %i' %
                  (part_counter, part_tot, usr_core_num))
        part_counter += 1

        # Start multiprocessing to get rank score
        # usr_core_num = 1
        if usr_core_num > 1:
            parallel_pool = Pool(usr_core_num)
            lipid_info_results_lst = []
            core_worker_count = 1
            for lipid_sub_lst in lipid_sub_key_lst:
                if isinstance(lipid_sub_lst, tuple) or isinstance(lipid_sub_lst, list):
                    if None in lipid_sub_lst:
                        lipid_sub_lst = filter(lambda x: x is not None, lipid_sub_lst)
                    else:
                        pass
                    if isinstance(lipid_sub_lst[0], tuple) or isinstance(lipid_sub_lst[0], list):
                        lipid_sub_dct = {k: lipid_spec_dct[k] for k in lipid_sub_lst}
                    else:
                        lipid_sub_dct = {lipid_sub_lst: lipid_spec_dct[lipid_sub_lst]}
                        lipid_sub_lst = tuple([lipid_sub_lst])
                    print('>>> >>> Core #%i ==> ...... processing ......' % core_worker_count)
                    if len(lipid_sub_dct.keys()) > 0:
                        lipid_info_result = parallel_pool.apply_async(get_lipid_info,
                                                                      args=(param_dct, usr_fa_df, checked_info_df,
                                                                            checked_info_groups, lipid_sub_lst,
                                                                            usr_weight_df, key_frag_dct,
                                                                            lipid_sub_dct, xic_dct))
                    lipid_info_results_lst.append(lipid_info_result)
                    core_worker_count += 1

            parallel_pool.close()
            parallel_pool.join()

            for lipid_info_result in lipid_info_results_lst:
                try:
                    tmp_lipid_info_df = lipid_info_result.get()
                except (KeyError, SystemError, ValueError):
                    tmp_lipid_info_df = 'error'
                    print('!!error!!--> This segment receive no Lipid identified.')
                if isinstance(tmp_lipid_info_df, str):
                    pass
                else:
                    if isinstance(tmp_lipid_info_df, pd.DataFrame):
                        if tmp_lipid_info_df.shape[0] > 0:
                            output_df = output_df.append(tmp_lipid_info_df)
        else:
            print('Using single core mode...')
            core_worker_count = 1
            for lipid_sub_lst in lipid_sub_key_lst:
                if isinstance(lipid_sub_lst, tuple) or isinstance(lipid_sub_lst, list):
                    if None in lipid_sub_lst:
                        lipid_sub_lst = filter(lambda x: x is not None, lipid_sub_lst)
                    else:
                        pass
                    if isinstance(lipid_sub_lst[0], tuple) or isinstance(lipid_sub_lst[0], list):
                        lipid_sub_dct = {k: lipid_spec_dct[k] for k in lipid_sub_lst}
                    else:
                        lipid_sub_dct = {lipid_sub_lst: lipid_spec_dct[lipid_sub_lst]}
                        lipid_sub_lst = tuple([lipid_sub_lst])
                    print('>>> Part %i Subset #%i ==> ...... processing ......' % (part_counter, core_worker_count))
                    if len(lipid_sub_dct.keys()) > 0:
                        tmp_lipid_info_df = get_lipid_info(param_dct, usr_fa_df, checked_info_df, checked_info_groups,
                                                           lipid_sub_lst, usr_weight_df, key_frag_dct,
                                                           lipid_sub_dct, xic_dct)

                    core_worker_count += 1
                    if isinstance(tmp_lipid_info_df, str):
                        pass
                    else:
                        if tmp_lipid_info_df.shape[0] > 0:
                            output_df = output_df.append(tmp_lipid_info_df)

    print('=== ==> --> Generate the output table')
    if output_df.shape[0] > 0:
        try:
            output_df = output_df.sort_values(by=['Lib_mz', 'Proposed_structures', 'MS2_scan_time', 'RANK_SCORE'])
        except KeyError:
            pass
        output_df.reset_index(drop=True, inplace=True)
        output_df.index += 1
        print('output_df')
        print(output_df.head(5))
        print(output_df.columns.tolist())
        print(output_df[['Proposed_structures', 'DISCRETE_ABBR', 'MS2_scan_time', 'img_name']])
        output_df.drop_duplicates(keep='first', inplace=True)
        log_pager.add_all_info(output_df)
        output_header_lst = output_df.columns.tolist()
        for _i_check in ['SN1_[FA-H]-_i', 'SN2_[FA-H]-_i', '[LPL(SN1)-H]-_i', '[LPL(SN2)-H]-_i',
                         '[LPL(SN1)-H2O-H]-_i', '[LPL(SN2)-H2O-H]-_i']:
            if _i_check not in output_header_lst:
                output_df[_i_check] = 0.0

        output_round_dct = {r'MS1_obs_mz': 4, r'Lib_mz': 4, 'ppm': 2, 'MS2_scan_time': 3,
                            'i_sn1': 2, 'i_sn2': 2, 'i_[M-H]-sn1': 2, 'i_[M-H]-sn2': 2,
                            'i_[M-H]-sn1-H2O': 2, 'i_[M-H]-sn2-H2O': 2
                            }
        # add intensities of target peaks to round list
        if len(target_ident_lst) > 0:
            for _t in target_ident_lst:
                output_round_dct[_t] = 2
        output_df = output_df.round(output_round_dct)

        output_df.rename(columns={'#Contaminated_peaks': '#Unspecific_peaks'}, inplace=True)

        output_short_lst = ['Proposed_structures', 'DISCRETE_ABBR', 'Formula_neutral', 'Formula_ion', 'Charge',
                            'Lib_mz', 'ppm', 'RANK_SCORE', 'MS1_obs_mz', 'MS1_obs_i', r'MS2_PR_mz', 'MS2_scan_time',
                            'DDA#', 'Scan#', 'SN1_[FA-H]-_i', 'SN2_[FA-H]-_i', '[LPL(SN1)-H]-_i', '[LPL(SN2)-H]-_i',
                            '[LPL(SN1)-H2O-H]-_i', '[LPL(SN2)-H2O-H]-_i']

        output_df = output_df[output_short_lst]
        output_df = output_df.sort_values(by=['MS1_obs_mz', 'MS2_scan_time', 'RANK_SCORE'],
                                          ascending=[True, True, False])
        output_df = output_df.reset_index(drop=True)
        output_df.index += 1
        try:
            output_df.to_excel(output_sum_xlsx, index=False)
            print(output_sum_xlsx)
        except IOError:
            output_df.to_excel('%s-%i%s' % (output_sum_xlsx[:-5], int(time.time()), '.xlsx'), index=False)
            print(output_sum_xlsx)
        print('=== ==> --> saved >>> >>> >>>')

    log_pager.close_page()

    tot_run_time = time.clock() - start_time

    # tot_run_time = '%.2f Sec\n' % tot_run_time

    print('>>> >>> >>> FINISHED in %s sec <<< <<< <<<' % tot_run_time)

    return tot_run_time, error_lst


if __name__ == '__main__':
    pl_class = 'PE'
    charge = '[M-H]-'
    # pl_class = 'PC'
    # charge = '[M+HCOO]-'
    mz_range = [650, 950]
    rt_range = [20, 30]
    count = 1

    usr_dct = {'fawhitelist_path_str': r'D:\project_lipidhunter\lipidhunterdev\ConfigurationFiles\FA_Whitelist.xlsx',
               'mzml_path_str': r'D:\project_lipidhunter\MF_mzML\MS2\070120_CM_neg_70min_SIN_I.mzML',
               'img_output_folder_str': r'D:\project_lipidhunter\lipidhunterdev\Temp\Test%s%i' % (pl_class, count),
               'xlsx_output_path_str': r'D:\project_lipidhunter\lipidhunterdev\Temp\Test%s%i\t%s_%i.xlsx'
                                       % (pl_class, count, pl_class, count),
               'lipid_specific_cfg':
                   r'D:\project_lipidhunter\lipidhunterdev\ConfigurationFiles\PL_specific_ion_cfg.xlsx',
               'hunter_start_time': '2017-12-21_15-27-49',
               'vendor': 'waters', 'experiment_mode': 'LC-MS', 'lipid_type': pl_class, 'charge_mode': charge,
               'rt_start': rt_range[0], 'rt_end': rt_range[1], 'mz_start': mz_range[0], 'mz_end': mz_range[1],
               'rank_score': True, 'rank_score_filter': 27.5, 'score_filter': 27.5,
               'isotope_score_filter': 75.0, 'fast_isotope': False,
               'ms_th': 1000, 'ms_ppm': 20, 'ms_max': 0, 'pr_window': 0.75, 'dda_top': 6,
               'ms2_th': 10, 'ms2_ppm': 50, 'ms2_infopeak_threshold': 0.001,
               'hg_th': 10.0, 'hg_ppm': 200.0, 'ms2_hginfopeak_threshold': 0.001,
               'score_cfg': r'D:\project_lipidhunter\lipidhunterdev\ConfigurationFiles\Score_cfg.xlsx',
               'fa_white_list_cfg': r'D:\project_lipidhunter\lipidhunterdev\ConfigurationFiles\FA_Whitelist.xlsx',
               'hunter_folder': r'D:\project_lipidhunter\lipidhunterdev',
               'core_number': 3, 'max_ram': 5, 'img_type': u'png', 'img_dpi': 300, 'tag_all_sn': True}
    log_lst = []
    t, log_lst = huntlipids(usr_dct, log_lst)
    print(t)
    print('test passed!')
