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
import multiprocessing
from multiprocessing import Pool
import os
import sys
from sys import platform
import time

from numpy import int64
import pandas as pd

try:
    from LibLipidHunter.LipidComposer import LipidComposer
    from LibLipidHunter.SpectraReader import extract_mzml
    from LibLipidHunter.SpectraReader import get_spectra
    from LibLipidHunter.SpectraReader import get_xic_from_pl
    from LibLipidHunter.SpectraReader import get_spec_info
    from LibLipidHunter.LogPageCreator import LogPageCreator
    from LibLipidHunter.PrecursorHunter import PrecursorHunter
    from LibLipidHunter.ScoreHunter import get_lipid_info
    from LibLipidHunter.PanelPlotter import gen_plot
except ImportError:  # for python 2.7.14
    from LipidComposer import LipidComposer
    from SpectraReader import extract_mzml
    from SpectraReader import get_spectra
    from SpectraReader import get_xic_from_pl
    from SpectraReader import get_spec_info
    from LogPageCreator import LogPageCreator
    from PrecursorHunter import PrecursorHunter
    from ScoreHunter import get_lipid_info
    from PanelPlotter import gen_plot


def huntlipids(param_dct, error_lst, save_fig=True):
    """

    :param(dict) param_dct:
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
                        'hunter_folder': r'D:\lipidhunter',
                        'core_number': 3, 'max_ram': 5, 'img_type': u'png', 'img_dpi': 300}
    :param(list) error_lst: empty list to store error info to display on GUI.
    :param(bool) save_fig: Can be set to False to skip image generation (not recommended).
    :return:
    """

    print('>>> hunter core start...')

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
    if 'xic_ppm' in list(param_dct.keys()):
        usr_xic_ppm = int(param_dct['xic_ppm'])
    else:
        usr_xic_ppm = usr_ms1_ppm

    # parameters from settings tab
    usr_core_num = param_dct['core_number']
    usr_max_ram = param_dct['max_ram']

    usr_dpi = param_dct['img_dpi']
    usr_img_type = param_dct['img_type']

    hunter_start_time_str = param_dct['hunter_start_time']

    if platform == "linux" or platform == "linux2":  # linux
        if usr_core_num > 1:
            os_typ = 'linux_multi'
            print('=== ==> --> LipidHunter Running on >>> Linux with multiprocessing mode ...')
        else:
            os_typ = 'linux_single'
            print('=== ==> --> LipidHunter Running on >>> Linux with single core mode ...')
    elif platform == "win32":  # Windows
        os_typ = 'windows'
        print('=== ==> --> LipidHunter Running on >>> Windows ...')
    else:
        usr_core_num = 1
        param_dct['core_number'] = 1
        os_typ = 'linux_single'
        print('=== ==> --> LipidHunter Running on >>> unoptimized system ... %s' % platform)
        print('!! Force to use single core mode !!')

    print('=== ==> --> Start to process >>>')
    print('=== ==> --> Lipid class: %s >>>' % usr_lipid_class)

    composer_param_dct = {'fa_whitelist': usr_fa_xlsx, 'lipid_type': usr_lipid_class,
                          'charge_mode': usr_charge, 'exact_position': 'FALSE'}

    existed_lipid_master_path = ''
    use_existed_lipid_master = False
    save_lipid_master_table = False
    if 'debug_mode' in list(param_dct.keys()):
        if param_dct['debug_mode'] == 'ON' and 'lipid_master_table' in list(param_dct.keys()):
            existed_lipid_master_path = param_dct['lipid_master_table']
            if os.path.isfile(existed_lipid_master_path):
                use_existed_lipid_master = True
            else:
                print('Failed to load existed Lipid Master table: %s', existed_lipid_master_path)

        if 'save_lipid_master_table' in list(param_dct.keys()):
            if param_dct['save_lipid_master_table'] == 'CSV':
                save_lipid_master_table = True

    if use_existed_lipid_master is False:
        try:
            print('==> --> Start to generate Lipid Master Table ...')
            usr_lipid_master_df = lipidcomposer.compose_lipid(param_dct=composer_param_dct, ms2_ppm=usr_ms2_ppm)
            print('=== ==> --> Lipid Master Table generated >>>', usr_lipid_master_df.shape[0])

        except Exception as e:
            print(e)
            error_lst.append(e)
            error_lst.append('Some files missing...')
            error_lst.append('Please check your settings in the configuration file ...')
            return False, error_lst, False
    else:
        try:
            print('Try to use existed Lipid Master table: %s' % existed_lipid_master_path)
            usr_lipid_master_df = pd.read_csv(existed_lipid_master_path)
            print('=== ==> --> Lipid Master table loaded >>>', usr_lipid_master_df.shape[0])
        except Exception as e:
            print(e)
            error_lst.append(e)
            return False, error_lst, False

    if save_lipid_master_table is True:
        log_master_name = 'Lipid_Master_%s.csv' % hunter_start_time_str
        log_master_name = os.path.join(output_folder, log_master_name)
        usr_lipid_master_df.to_csv(log_master_name)
        print('==> --> Lipid Master table Saved as: ', log_master_name)
    else:
        pass

    # for TG has the fragment of neutral loss of the FA and the fragments for the MG
    usr_fa_df = lipidcomposer.calc_fa_query(usr_lipid_class, usr_fa_xlsx, ms2_ppm=usr_ms2_ppm)

    lipid_info_df = usr_lipid_master_df

    # cut lib info to the user defined m/z range
    # TODO (georgia.angelidou@uni-leipzig.de): support for the sphingomyelins and ceramides
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
            if usr_charge == '[M+NH4]+':
                lipid_info_df = lipid_info_df[
                    (mz_start <= lipid_info_df['[M+NH4]+_MZ']) & (lipid_info_df['[M+NH4]+_MZ'] <= mz_end)]
            elif usr_charge == '[M+H]+':
                lipid_info_df = lipid_info_df[
                    (mz_start <= lipid_info_df['[M+H]+_MZ']) & (lipid_info_df['[M+H]+_MZ'] <= mz_end)]
            elif usr_charge == '[M+Na]+':
                lipid_info_df = lipid_info_df[
                    (mz_start <= lipid_info_df['[M+Na]+_MZ']) & (lipid_info_df['[M+Na]+_MZ'] <= mz_end)]
    else:
        error_lst.append('Lipid class or charge NOT supported.  User input lipid class = %s, charge = %s. '
                         % (usr_lipid_class, usr_charge))

    print('=== ==> --> Lipid Master table Ready ...')

    # TODO(zhixu.ni@uni-leipzig.de): Add more error to the error_lst.

    pr_hunter = PrecursorHunter(lipid_info_df, param_dct, os_type=os_typ)

    output_df = pd.DataFrame()

    print('=== ==> --> Start to process')
    print('=== ==> --> Lipid class: %s' % usr_lipid_class)

    # generate the Weight factor df
    usr_weight_df = pd.read_excel(score_cfg, index_col='Type')

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

    if ms1_obs_pr_df is False:
        print('!! NO suitable precursor --> Check settings!!\n')
        error_lst.append('!! NO suitable precursor --> Check settings!!\n')
        return False, error_lst, False

    print('=== ==> --> ms1 precursor matched')

    # remove bad precursors, keep the matched scans by DDA_rank and scan number
    # build unique identifier for each scan with scan_number00dda_rank use numpy.int64 to avoid large scan_number
    # usr_scan_info_df['scan_checker'] = ((10000 * usr_scan_info_df['scan_number'] + usr_scan_info_df['DDA_rank'])
    #                                     .astype(int64))
    # ms1_obs_pr_df['scan_checker'] = ((10000 * ms1_obs_pr_df['scan_number'] + ms1_obs_pr_df['DDA_rank'])
    #                                  .astype(int64))
    usr_scan_info_df['scan_checker'] = (usr_scan_info_df['scan_number'].astype(int64).astype(str).str.
                                        cat(usr_scan_info_df['DDA_rank'].astype(int64).astype(str), sep='_'))
    ms1_obs_pr_df['scan_checker'] = (ms1_obs_pr_df['scan_number'].astype(int64).astype(str).str.
                                     cat(ms1_obs_pr_df['DDA_rank'].astype(int64).astype(str).astype(str), sep='_'))

    usr_scan_checker_lst = usr_scan_info_df['scan_checker'].tolist()
    checked_info_df = ms1_obs_pr_df[ms1_obs_pr_df['scan_checker'].isin(usr_scan_checker_lst)]
    checked_info_df.is_copy = False
    checked_info_df.sort_values(by=['scan_checker', 'Lib_mz'], ascending=[True, True], inplace=True)
    if 'debug_mode' in list(param_dct.keys()):
        if param_dct['debug_mode'] == 'ON':
            usr_scan_info_df.to_csv(os.path.join(output_folder, 'usr_scan_info.csv'))
            ms1_obs_pr_df.to_csv(os.path.join(output_folder, 'ms1_obs_pr_df.csv'))
            checked_info_df.to_csv(os.path.join(output_folder, 'checked_info_df.csv'))

    if checked_info_df.shape[0] == 0:
        print('!! No identification in pre-match steps !!')
        error_lst.append('!! No identification in pre-match steps !!\n')
        return False, error_lst, False
    else:
        print('>>> features identified in the pre-match: ', checked_info_df.shape[0])

    ms1_xic_mz_lst = ms1_obs_pr_df['MS1_XIC_mz'].values.tolist()
    ms1_xic_mz_lst = sorted(set(ms1_xic_mz_lst))
    print('ms1_xic_mz_lst', len(ms1_xic_mz_lst))
    print(ms1_xic_mz_lst)
    print('=== ==> --> Start to extract XIC')
    if len(ms1_xic_mz_lst) >= 3 * usr_core_num:
        sub_len = int(math.ceil(len(ms1_xic_mz_lst) / usr_core_num))
        # core_key_list = list(*(iter(ms1_xic_mz_lst),) * sub_len)
        core_key_list = [ms1_xic_mz_lst[k: k + sub_len] for k in range(0, len(ms1_xic_mz_lst), sub_len)]

    else:
        core_key_list = [ms1_xic_mz_lst]

    # Start multiprocessing to get XIC
    print('!!!!!! Start multiprocessing to get XIC ==> ==> ==> Number of Cores: %i' % usr_core_num)
    xic_dct = {}

    if usr_core_num > 1:
        xic_results_lst = []

        if os_typ == 'windows':
            parallel_pool = Pool(usr_core_num)
            queue = ''
            worker_count = 1
            for core_list in core_key_list:
                if isinstance(core_list, tuple) or isinstance(core_list, list):
                    if None in core_list:
                        core_list = [x for x in core_list if x is not None]
                    else:
                        pass
                    print('>>> >>> Core #%i ==> ...... processing ......' % worker_count)
                    print(core_list)
                    xic_result = parallel_pool.apply_async(get_xic_from_pl, args=(core_list, ms1_xic_df, usr_xic_ppm,
                                                                                  os_typ, queue))
                    worker_count += 1
                    xic_results_lst.append(xic_result)

            parallel_pool.close()
            parallel_pool.join()

            for xic_result in xic_results_lst:
                try:
                    sub_xic_dct = xic_result.get()
                    if len(list(sub_xic_dct.keys())) > 0:
                        xic_dct.update(sub_xic_dct)

                except (KeyError, SystemError, ValueError):
                    pass

        else:  # for linux
            jobs = []
            queue = multiprocessing.Queue()
            worker_count = 1
            for core_list in core_key_list:
                if isinstance(core_list, tuple) or isinstance(core_list, list):
                    if None in core_list:
                        core_list = [x for x in core_list if x is not None]
                    else:
                        pass
                    print('>>> >>> Core #%i ==> ...... processing ......' % worker_count)
                    print(core_list)
                    job = multiprocessing.Process(target=get_xic_from_pl, args=(core_list, ms1_xic_df, usr_xic_ppm,
                                                                                os_typ, queue))
                    worker_count += 1
                    jobs.append(job)
                    job.start()
                    xic_results_lst.append(queue.get())
            for j in jobs:
                j.join()

            for xic_result in xic_results_lst:
                try:
                    if len(list(xic_result.keys())) > 0:
                        xic_dct.update(xic_result)
                except (KeyError, SystemError, ValueError):
                    pass
    else:
        print('Using single core mode...')
        queue = ''
        worker_count = 1
        for core_list in core_key_list:
            if isinstance(core_list, tuple) or isinstance(core_list, list):
                if None in core_list:
                    core_list = [x for x in core_list if x is not None]
                else:
                    pass
                print('>>> >>> Core #%i ==> ...... processing ......' % worker_count)
                print(core_list)
                sub_xic_dct = get_xic_from_pl(core_list, ms1_xic_df, usr_xic_ppm, os_typ, queue)
                worker_count += 1
                if len(list(sub_xic_dct.keys())) > 0:
                    xic_dct.update(sub_xic_dct)

    # print('xic_dct', len(xic_dct.keys()))
    # print(xic_dct.keys())

    if len(list(xic_dct.keys())) == 0:
        print('!! No precursor for XIC found !!')
        error_lst.append('!! No precursor for XIC found !!\n')
        return False, error_lst, False
    else:
        print('=== ==> --> Number of XIC extracted: %i' % len(list(xic_dct.keys())))

    target_ident_lst = []

    print('=== ==> --> Start to Hunt for Lipids !!')
    # checked_info_groups = checked_info_df.groupby(['Lib_mz', 'MS2_PR_mz', 'Formula', 'scan_time', 'Ion'])
    checked_info_groups = checked_info_df.groupby(['Formula', 'scan_checker'])
    lipid_all_group_key_lst = list(checked_info_groups.groups.keys())
    spec_sub_len = int(math.ceil(len(lipid_all_group_key_lst) / usr_core_num))
    spec_sub_key_lst = [lipid_all_group_key_lst[k: k + spec_sub_len] for k in range(0, len(lipid_all_group_key_lst),
                                                                                    spec_sub_len)]
    lipid_spec_info_dct = {}
    # TODO (georiga.angelidou@uni-leipzig.de): remember to remove
    # spec_sub_key_lst = [[(847.67916, 847.680061, 'C53H92O6Na+', 22.063973, '[M+Na]+'),
    #                      (847.67916, 847.679436, 'C53H92O6Na+', 22.598964, '[M+Na]+'),
    #                      (847.67916, 847.679513, 'C53H92O6Na+', 22.409925, '[M+Na]+'),
    #                      (847.67916, 847.67998, 'C53H92O6Na+', 22.235679, '[M+Na]+')]]
    if usr_core_num > 1:
        spec_results_lst = []

        if os_typ == 'windows':
            parallel_pool = Pool(usr_core_num)
            queue = ''
            worker_count = 1
            for _sub_lst in spec_sub_key_lst:
                if isinstance(_sub_lst, tuple) or isinstance(_sub_lst, list):
                    if None in _sub_lst:
                        _sub_lst = [x for x in _sub_lst if x is not None]
                    else:
                        pass
                    print('>>> >>> Core #%i ==> ...... processing ......' % worker_count)
                    spec_result = parallel_pool.apply_async(get_spec_info, args=(_sub_lst, checked_info_groups,
                                                                                 usr_scan_info_df, os_typ, queue))
                    worker_count += 1
                    spec_results_lst.append(spec_result)

            parallel_pool.close()
            parallel_pool.join()

            for spec_result in spec_results_lst:
                try:
                    sub_spec_dct = spec_result.get()
                    print(sub_spec_dct)
                    if len(list(sub_spec_dct.keys())) > 0:
                        lipid_spec_info_dct.update(sub_spec_dct)
                except (KeyError, SystemError, ValueError):
                    print('ValueError: must supply a tuple to get_group with multiple grouping keys')

        else:  # for linux
            jobs = []
            queue_spec = multiprocessing.Queue()
            worker_count = 1
            for _sub_lst in spec_sub_key_lst:
                if isinstance(_sub_lst, tuple) or isinstance(_sub_lst, list):
                    if None in _sub_lst:
                        _sub_lst = [x for x in _sub_lst if x is not None]
                    else:
                        pass
                    print('>>> >>> Core #%i ==> ...... processing ......' % worker_count)
                    job = multiprocessing.Process(target=get_spec_info, args=(_sub_lst, checked_info_groups,
                                                                              usr_scan_info_df, os_typ, queue_spec))
                    worker_count += 1
                    jobs.append(job)
                    job.start()
                    spec_results_lst.append(queue_spec.get())
            for j in jobs:
                j.join()

            for spec_result in spec_results_lst:
                try:
                    if len(list(spec_result.keys())) > 0:
                        lipid_spec_info_dct.update(spec_result)
                except (KeyError, SystemError, ValueError):
                    print('ValueError: must supply a tuple to get_group with multiple grouping keys')
    else:
        print('Using single core mode...')
        queue = ''
        worker_count = 1
        for _sub_lst in spec_sub_key_lst:
            if isinstance(_sub_lst, tuple) or isinstance(_sub_lst, list):
                if None in _sub_lst:
                    _sub_lst = [x for x in _sub_lst if x is not None]
                else:
                    if isinstance(_sub_lst[0], float):
                        _sub_lst3 = ()
                        _sub_lst3 = _sub_lst3 + (_sub_lst,)
                        _sub_lst = _sub_lst3

                print('>>> >>> Core #%i ==> ...... processing ......' % worker_count)
                sub_spec_dct = get_spec_info(_sub_lst, checked_info_groups, usr_scan_info_df, os_typ, queue)
                worker_count += 1
                if len(list(sub_spec_dct.keys())) > 0:
                    lipid_spec_info_dct.update(sub_spec_dct)

            else:
                pass

    print('lipid_spec_info_dct', len(list(lipid_spec_info_dct.keys())))

    # Single process ONLY. usr_spectra_pl is too big in RAM --> RAM leaking during copy
    lipid_spec_dct = {}
    spec_info_key_lst = list(lipid_spec_info_dct.keys())
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

    found_spec_key_lst = list(lipid_spec_dct.keys())
    found_spec_key_lst = sorted(found_spec_key_lst, key=lambda x: x[0])
    spec_key_num = len(found_spec_key_lst)
    print('spec_key_num', spec_key_num)

    # parse specific peak info
    pl_class_lst = ['PA', 'PC', 'PE', 'PG', 'PI', 'PS', 'PIP']
    pl_neg_chg_lst = ['[M-H]-', '[M+HCOO]-', '[M+CH3COO]-']
    tg_class_lst = ['TG', 'DG']
    tg_pos_chg_lst = ['[M+NH4]+', '[M+H]+', '[M+Na]+']
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
    elif usr_lipid_class in tg_class_lst and usr_charge in tg_pos_chg_lst:
        charge_mode = 'POS'
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

    # Start multiprocessing to get rank score

    lipid_info_img_lst = []

    if usr_core_num > 1:

        lipid_part_key_lst = []
        split_seg = 1

        if spec_key_num > (usr_core_num * 24):

            # Split tasks into few parts to avoid core waiting in multiprocessing
            if usr_core_num * 24 < spec_key_num <= usr_core_num * 48:
                split_seg = 2
            elif usr_core_num * 48 < spec_key_num <= usr_core_num * 96:
                split_seg = 3
            elif usr_core_num * 96 < spec_key_num:
                split_seg = 4
            else:
                split_seg = 1

            lipid_part_len = int(math.ceil(spec_key_num / split_seg))
            lipid_part_lst = [found_spec_key_lst[k: k + lipid_part_len] for k in range(0, spec_key_num,
                                                                                       lipid_part_len)]
            print('lipid_part_number: ', len(lipid_part_lst), ' lipid_part_len:', lipid_part_len)

            for part_lst in lipid_part_lst:
                if None in part_lst:
                    part_lst = [x for x in part_lst if x is not None]
                lipid_sub_len = int(math.ceil(len(part_lst) / usr_core_num))
                # print('lipid_sub_len', lipid_sub_len)
                lipid_sub_key_lst = [part_lst[k: k + lipid_sub_len] for k in range(0, len(part_lst), lipid_sub_len)]
                lipid_part_key_lst.append(lipid_sub_key_lst)

        else:
            lipid_sub_len = int(math.ceil(spec_key_num / usr_core_num))
            lipid_sub_key_lst = [found_spec_key_lst[k: k + lipid_sub_len] for k in
                                 range(0, spec_key_num, lipid_sub_len)]
            lipid_part_key_lst.append(lipid_sub_key_lst)

        # print('lipid_part_number: ', len(lipid_part_key_lst), ' lipid_part_len:', len(lipid_part_key_lst[0]))

        part_tot = len(lipid_part_key_lst)
        # print('part_tot', part_tot)
        # print(lipid_part_key_lst)
        part_counter = 1
        queue = ''

        for lipid_sub_key_lst in lipid_part_key_lst:

            lipid_info_results_lst = []

            if part_tot == 1:
                print('>>> Start multiprocessing to get Score ==> Max Number of Cores: %i' % usr_core_num)
            else:
                print('>>> Start multiprocessing to get Score ==> Part %i / %i '
                      '--> Max Number of Cores: %i | x%i Features each'
                      % (part_counter, part_tot, usr_core_num, split_seg))

            if os_typ == 'windows':
                parallel_pool = Pool(usr_core_num)

                worker_count = 1
                for lipid_sub_lst in lipid_sub_key_lst:
                    if isinstance(lipid_sub_lst, tuple) or isinstance(lipid_sub_lst, list):
                        if None in lipid_sub_lst:
                            lipid_sub_lst = [x for x in lipid_sub_lst if x is not None]
                        else:
                            pass
                        if isinstance(lipid_sub_lst[0], tuple) or isinstance(lipid_sub_lst[0], list):
                            lipid_sub_dct = {k: lipid_spec_dct[k] for k in lipid_sub_lst}
                        else:
                            lipid_sub_dct = {lipid_sub_lst: lipid_spec_dct[lipid_sub_lst]}
                            lipid_sub_lst = tuple([lipid_sub_lst])
                        print('>>> >>> Core #%i ==> ...... processing ......' % worker_count)
                        if len(list(lipid_sub_dct.keys())) > 0:
                            lipid_info_result = parallel_pool.apply_async(get_lipid_info,
                                                                          args=(param_dct, usr_fa_df, checked_info_df,
                                                                                checked_info_groups, lipid_sub_lst,
                                                                                usr_weight_df, key_frag_dct,
                                                                                lipid_sub_dct, xic_dct, worker_count,
                                                                                save_fig, os_typ, queue))
                            lipid_info_results_lst.append(lipid_info_result)
                            worker_count += 1

                parallel_pool.close()
                parallel_pool.join()

            else:  # for linux
                jobs = []
                queue = multiprocessing.Queue()
                worker_count = 1
                for lipid_sub_lst in lipid_sub_key_lst:
                    if isinstance(lipid_sub_lst, tuple) or isinstance(lipid_sub_lst, list):
                        if None in lipid_sub_lst:
                            lipid_sub_lst = [x for x in lipid_sub_lst if x is not None]
                        else:
                            pass
                        if isinstance(lipid_sub_lst[0], tuple) or isinstance(lipid_sub_lst[0], list):
                            lipid_sub_dct = {k: lipid_spec_dct[k] for k in lipid_sub_lst}
                        else:
                            lipid_sub_dct = {lipid_sub_lst: lipid_spec_dct[lipid_sub_lst]}
                            lipid_sub_lst = tuple([lipid_sub_lst])
                        print('>>> >>> Core #%i ==> ...... processing ......' % worker_count)
                        job = multiprocessing.Process(target=get_lipid_info, args=(param_dct, usr_fa_df,
                                                                                   checked_info_df,
                                                                                   checked_info_groups, lipid_sub_lst,
                                                                                   usr_weight_df, key_frag_dct,
                                                                                   lipid_sub_dct, xic_dct, worker_count,
                                                                                   save_fig, os_typ, queue))
                        worker_count += 1
                        jobs.append(job)
                        job.start()
                        lipid_info_results_lst.append(queue.get())
                for j in jobs:
                    j.join()

            # Merge multiprocessing results
            for lipid_info_result in lipid_info_results_lst:

                if os_typ == 'windows':
                    try:
                        tmp_lipid_info = lipid_info_result.get()
                        tmp_lipid_info_df = tmp_lipid_info[0]
                        # print(tmp_lipid_info_df)
                        tmp_lipid_img_lst = tmp_lipid_info[1]
                        # print('tmp_lipid_img_lst')
                        # print(len(tmp_lipid_img_lst))
                    except (KeyError, SystemError, ValueError, TypeError):
                        tmp_lipid_info_df = 'error'
                        tmp_lipid_img_lst = []
                        print('!!error!!--> This segment receive no Lipid identified.')
                else:  # for linux
                    try:
                        tmp_lipid_info_df = lipid_info_result[0]
                        tmp_lipid_img_lst = lipid_info_result[1]

                    except (KeyError, SystemError, ValueError, TypeError):
                        tmp_lipid_info_df = 'error'
                        tmp_lipid_img_lst = []
                        print('!!error!!--> This segment receive no Lipid identified.')

                if isinstance(tmp_lipid_info_df, str):
                    pass
                else:
                    if isinstance(tmp_lipid_info_df, pd.DataFrame):
                        if tmp_lipid_info_df.shape[0] > 0:
                            output_df = output_df.append(tmp_lipid_info_df)
                            lipid_info_img_lst.extend(tmp_lipid_img_lst)
            if part_tot == 1:
                print('>>> multiprocessing results merged')
            else:
                print('>>> multiprocessing results merged ==> Part %i / %i '
                      % (part_counter, part_tot))
            part_counter += 1

    else:
        print('Using single core mode...')
        worker_count = 1
        lipid_info_results_lst = get_lipid_info(param_dct, usr_fa_df, checked_info_df,
                                                checked_info_groups, found_spec_key_lst, usr_weight_df,
                                                key_frag_dct, lipid_spec_dct, xic_dct, worker_count)
        tmp_lipid_info_df = lipid_info_results_lst[0]
        tmp_lipid_img_lst = lipid_info_results_lst[1]
        if isinstance(tmp_lipid_info_df, pd.DataFrame):
            if tmp_lipid_info_df.shape[0] > 0:
                output_df = output_df.append(tmp_lipid_info_df)
                lipid_info_img_lst = tmp_lipid_img_lst

    print('=== ==> --> Generate the output table')
    if output_df.shape[0] > 0:
        try:
            output_df = output_df.sort_values(by=['Lib_mz', 'Bulk_identification', 'MS2_scan_time', 'RANK_SCORE'],
                                              ascending=[True, True, True, False])
        except KeyError:
            pass
        output_df.reset_index(drop=True, inplace=True)
        output_df.index += 1
        # print('output_df')
        # print(output_df.head(5))
        # print(output_df.columns.values.tolist())
        # print(output_df[['Proposed_structures', 'DISCRETE_ABBR', 'MS2_scan_time', 'img_name']])

        output_df.drop_duplicates(keep='first', inplace=True)
        output_header_lst = output_df.columns.values.tolist()
        # TODO (georgia.angeldou@uni-leipzig.de): Add the info for the DG
        # TODO (zhixu.ni@uni-leipzig.de): check following if segment
        if usr_lipid_class in ['PA', 'PC', 'PE', 'PG', 'PI', 'PIP', 'PS']:
            output_list = ['FA1_[FA-H]-_i', 'FA2_[FA-H]-_i', '[LPL(FA1)-H]-_i', '[LPL(FA2)-H]-_i',
                           '[LPL(FA1)-H2O-H]-_i', '[LPL(FA2)-H2O-H]-_i']
            output_round_dct = {r'MS1_obs_mz': 4, r'Lib_mz': 4, 'ppm': 2, 'MS2_scan_time': 3,
                                'i_fa1': 2, 'i_fa2': 2, 'i_[M-H]-fa1': 2, 'i_[M-H]-fa2': 2,
                                'i_[M-H]-fa1-H2O': 2, 'i_[M-H]-fa2-H2O': 2
                                }
            # for _i_check in ['SN1_[FA-H]-_i', 'SN2_[FA-H]-_i', '[LPL(SN1)-H]-_i', '[LPL(SN2)-H]-_i',
            #                  '[LPL(SN1)-H2O-H]-_i', '[LPL(SN2)-H2O-H]-_i']:
            #     if _i_check not in output_header_lst:
            #         output_df[_i_check] = 0.0

        elif usr_lipid_class in ['TG'] and usr_charge in ['[M+NH4]+', '[M+H]+']:
            output_list = ['FA1_[FA-H2O+H]+_i', 'FA2_[FA-H2O+H]+_i', 'FA3_[FA-H2O+H]+_i', '[MG(FA1)-H2O+H]+_i',
                           '[MG(FA2)-H2O+H]+_i', '[MG(FA3)-H2O+H]+_i', '[M-(FA1)+H]+', '[M-(FA2)+H]+_i',
                           '[M-(FA3)+H]+_i']
            output_round_dct = {r'MS1_obs_mz': 4, r'Lib_mz': 4, 'ppm': 2, 'MS2_scan_time': 3,
                                'i_fa1': 2, 'i_fa2': 2, 'i_fa3': 2, 'i_[M+H]-fa1': 2, 'i_[M+H]-fa2': 2,
                                'i_[M+H]-fa3': 2, 'i_[MG(fa1)+H]-H2O': 2, 'i_[MG(fa2)+H]-H2O': 2,
                                'i_[MG(fa3)+H]-H2O': 2}
        elif usr_lipid_class in ['TG'] and usr_charge in ['[M+Na]+']:
            # TODO (georgia.angelidou@uni-leipzig.de): check why this can cause some problems with [MGSN1-H2O+H]+
            output_list = ['FA1_[FA-H2O+H]+_i', 'FA2_[FA-H2O+H]+_i', 'FA3_[FA-H2O+H]+_i', '[MG(FA1)-H2O+H]+_i',
                           '[MG(FA3)-H2O+H]+_i', '[MG(FA3)-H2O+H]+_i', '[M-(FA1)+Na]+_i', '[M-(FA2)+Na]+_i',
                           '[M-(FA3)+Na]+_i', '[M-(FA1-H+Na)+N]+_i', '[M-(FA2-H+Na)+H]+_i']
            output_round_dct = {r'MS1_obs_mz': 4, r'lib_mz': 4, 'ppm': 2, 'MS2_scan_time': 3, 'i_fa1': 2, 'i_fa2': 2,
                                'i_fa3': 2, 'i_[M+Na]-fa1': 2, 'i_[M+Na]-fa2': 2, 'i_[M+Na]-fa3': 2,
                                'i_[M+H]-fa1-H+Na': 2, 'i_[M+H]-fa2-H+Na': 2, 'i_[M+H]-fa3-H+Na': 2}
        elif usr_lipid_class in ['DG'] and usr_charge in ['[M+H]+', '[M+NH4]+', '[M+Na]+']:
            output_list = ['FA1_[FA-H2O+H]+_i', 'FA2_[FA-H2O+H]+_i', '[MG(FA1)-H2O+H]+_i',
                           '[MG(FA2)-H2O+H]+_i']
            output_round_dct = {r'MS1_obs_mz': 4, r'Lib_mz': 4, 'ppm': 2, 'MS2_scan_time': 3, 'i_fa1': 2, 'i_fa2': 2,
                                'i_[MG(fa1)+H]-H2O': 2, 'i_[MG(fa2)+H]-H2O': 2}

        else:
            output_list = ['FA1_[FA-H]-_i', 'FA2_[FA-H]-_i', '[LPL(FA1)-H]-_i', '[LPL(FA2)-H]-_i',
                           '[LPL(FA1)-H2O-H]-_i', '[LPL(FA2)-H2O-H]-_i']
            output_round_dct = {r'MS1_obs_mz': 4, r'Lib_mz': 4, 'ppm': 2, 'MS2_scan_time': 3,
                                'i_fa1': 2, 'i_fa2': 2, 'i_[M-H]-fa1': 2, 'i_[M-H]-fa2': 2,
                                'i_[M-H]-fa1-H2O': 2, 'i_[M-H]-fa2-H2O': 2
                                }
        for _i_check in output_list:
            if _i_check not in output_header_lst:
                output_df[_i_check] = 0.0

        # add intensities of target peaks to round list
        if len(target_ident_lst) > 0:
            for _t in target_ident_lst:
                output_round_dct[_t] = 2
        output_df = output_df.round(output_round_dct)

        output_df.rename(columns={'#Contaminated_peaks': '#Unspecific_peaks'}, inplace=True)
        # TODO (georgia.angelidou@uni-leipzig.de): Add the DG section
        if usr_lipid_class in ['TG'] and usr_charge in ['[M+H]+', '[M+NH4]+']:
            output_short_lst = ['Proposed_structures', 'DISCRETE_ABBR', 'Formula_neutral', 'Formula_ion', 'Charge',
                                'Lib_mz', 'ppm', 'RANK_SCORE', 'MS1_obs_mz', 'MS1_obs_i', r'MS2_PR_mz', 'MS2_scan_time',
                                'DDA#', 'Scan#', 'FA1_[FA-H2O+H]+_i', 'FA2_[FA-H2O+H]+_i', 'FA3_[FA-H2O+H]+_i',
                                '[MG(FA1)-H2O+H]+_i', '[MG(FA2)-H2O+H]+_i', '[MG(FA3)-H2O+H]+_i', '[M-(FA1)+H]+_i',
                                '[M-(FA2)+H]+_i', '[M-(FA3)+H]+_i']
        elif usr_lipid_class in ['TG'] and usr_charge in ['[M+Na]+']:
            # TODO (georgia.angelidou@uni-leipzig.de): need to solve the problem for the below sections
            # '[MG(SN1)-H2O+H]+_i', '[MG(SN2)-H2O+H]+_i', '[MG(SN3)-H2O+H]+_i',
            output_short_lst = ['Proposed_structures', 'DISCRETE_ABBR', 'Formula_neutral', 'Formula_ion', 'Charge',
                                'Lib_mz', 'ppm', 'RANK_SCORE', 'MS1_obs_mz', 'MS1_obs_i', r'MS2_PR_mz', 'MS2_scan_time',
                                'DDA#', 'Scan#', 'FA1_[FA-H2O+H]+_i', 'FA2_[FA-H2O+H]+_i', 'FA3_[FA-H2O+H]+_i',
                                '[M-(FA1)+Na]+_i',
                                '[M-(FA2)+Na]+_i', '[M-(FA3)+Na]+_i']
        elif usr_lipid_class in ['DG'] and usr_charge in ['[M+H]+', '[M+NH4]+', '[M+Na]+']:
            # problem with the following key:  'SN2_[FA-H2O+H]_i',
            output_short_lst = ['Proposed_structures', 'DISCRETE_ABBR', 'Formula_neutral', 'Formula_ion', 'Charge',
                                'Lib_mz', 'ppm', 'RANK_SCORE', 'MS1_obs_mz', 'MS1_obs_i', r'MS2_PR_mz', 'MS2_scan_time',
                                'DDA#', 'Scan#', 'FA1_[FA-H2O+H]+_i', '[MG(FA1)-H2O+H]+_i',
                                '[MG(FA2)-H2O+H]+_i']
        else:
            output_short_lst = ['Proposed_structures', 'DISCRETE_ABBR', 'Formula_neutral', 'Formula_ion', 'Charge',
                                'Lib_mz', 'ppm', 'RANK_SCORE', 'MS1_obs_mz', 'MS1_obs_i', r'MS2_PR_mz', 'MS2_scan_time',
                                'DDA#', 'Scan#', 'FA1_[FA-H]-_i', 'FA2_[FA-H]-_i', '[LPL(FA1)-H]-_i', '[LPL(FA2)-H]-_i',
                                '[LPL(FA1)-H2O-H]-_i', '[LPL(FA2)-H2O-H]-_i']

        final_output_df = output_df[output_short_lst]
        final_output_df = final_output_df.sort_values(by=['MS1_obs_mz', 'MS2_scan_time', 'RANK_SCORE'],
                                                      ascending=[True, True, False])
        final_output_df = final_output_df.reset_index(drop=True)
        final_output_df.index += 1

        output_sum_xlsx_directory = os.path.dirname(output_sum_xlsx)
        if not os.path.exists(output_sum_xlsx_directory):
            os.makedirs(output_sum_xlsx_directory)
        try:
            final_output_df.to_excel(output_sum_xlsx, index=False)
            print(output_sum_xlsx)
        except IOError:
            final_output_df.to_excel('%s-%i%s' % (output_sum_xlsx[:-5], int(time.time()), '.xlsx'), index=False)
            print(output_sum_xlsx)
        print('=== ==> --> saved >>> >>> >>>')

    else:
        error_lst.append('!! Warning !! NO Lipid identified in this file.\n!! Please check your settings !!')
        tot_run_time = time.clock() - start_time
        print('!! Warning !!--> This file got no Lipid identified.')
        print('=== == --> identification finished in %s sec <<< <<< <<<' % tot_run_time)
        return tot_run_time, error_lst, output_df

    # Start multiprocessing to save img
    if save_fig is True:

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

        # generate html files
        log_pager = LogPageCreator(output_folder, hunter_start_time_str, param_dct)
        log_pager.add_all_info(output_df)
        log_pager.close_page()

        print('=== == --> start to generate images: image count %i' % len(lipid_info_img_lst))

        if usr_core_num > 1:
            parallel_pool = Pool(usr_core_num)
            img_num = len(lipid_info_img_lst)
            img_sub_len = int(math.ceil(img_num / usr_core_num))
            img_sub_key_lst = [lipid_info_img_lst[k: k + img_sub_len] for k in range(0, img_num, img_sub_len)]

            worker_count = 1
            for img_sub_lst in img_sub_key_lst:
                if isinstance(img_sub_lst, tuple) or isinstance(img_sub_lst, list):
                    if None in img_sub_lst:
                        img_sub_lst = [x for x in img_sub_lst if x is not None]
                    else:
                        pass
                    # img_params_dct = {'lipid_info_img_lst': img_sub_lst, 'usr_core_num': usr_core_num,
                    #                   'usr_img_type': usr_img_type, 'usr_dpi': usr_dpi, 'usr_vendor': usr_vendor,
                    #                   'usr_ms1_precision': usr_ms1_precision, 'worker_count': worker_count}

                    if len(img_sub_lst) > 0:
                        print('>>> >>> Core #%i ==> Generating output images ... image count: %i'
                              % (worker_count, len(img_sub_lst)))
                        if 'debug_mode' in list(param_dct.keys()):
                            if param_dct['debug_mode'] == 'ON':
                                for img_param_dct in img_sub_lst:
                                    print(img_param_dct['save_img_as'])
                        parallel_pool.apply_async(gen_plot, args=(img_sub_lst, worker_count, usr_img_type,
                                                                  usr_dpi, usr_vendor, usr_ms1_precision))
                        worker_count += 1

            parallel_pool.close()
            parallel_pool.join()

        else:
            print('Using single core mode...')
            if isinstance(lipid_info_img_lst, tuple) or isinstance(lipid_info_img_lst, list):
                if None in lipid_info_img_lst:
                    lipid_info_img_lst = [x for x in lipid_info_img_lst if x is not None]
                else:
                    pass
                if len(lipid_info_img_lst) > 0:
                    gen_plot(lipid_info_img_lst, worker_count, usr_img_type, usr_dpi,
                             usr_vendor, usr_ms1_precision)
    else:
        print('!!!!!! User skip image generation !!!!!!')

    tot_run_time = time.clock() - start_time

    print('>>> >>> >>> FINISHED in %s sec <<< <<< <<<' % tot_run_time)

    return tot_run_time, error_lst, output_df


if __name__ == '__main__':

    # set the core number and max ram in GB to be used for the test
    core_count = 1
    max_ram = 5  # int only
    save_images = True  # True --> generate images, False --> NO images (not recommended)

    # full_test_lst = [['PC', 'waters'],['PE', 'waters'], ['TG', 'waters','[M+H]+'], ['TG', 'waters', '[M+NH4]+'],
    # ['TG', 'thermo', '[M+NH4]+']]

    usr_test_lst = [
        ['TG', 'thermo', '[M+NH4]+', 'TG_waters_NH4'],
        # ['PC', 'waters', '[M+HCOO]-', 'PC_waters'],
    ]

    # set the default files
    pl_mzml_waters = r'../Test/mzML/PL_neg_waters_synapt-g2si.mzML'  # Ni file
    tg_mzml_waters = r'../Test/mzML/TG_pos_waters_synapt-g2si.mzML'  # Mile file
    tg_mzml_thermo = r'../Test/mzML/TG_pos_thermo_2.mzML'  # Angela
    tg_mzml_SCIEXS = r'../Test/mzML/20140613_HSL002_Positive_01.mzML'  # Dataset
    tg_mzml_agilent = r'../Test/mzML/Test_agilent.mzML'  # position holder

    pl_base_dct = {'fawhitelist_path_str': r'../ConfigurationFiles/01-FA_Whitelist_PL.xlsx',
                   'lipid_specific_cfg': r'../ConfigurationFiles/02-Specific_ions.xlsx',
                   'score_cfg': r'../ConfigurationFiles/03-Score_weight_PL.xlsx'}

    tg_base_dct = {'fawhitelist_path_str': r'../ConfigurationFiles/01-FA_Whitelist_TG-DG.xlsx',
                   'lipid_specific_cfg': r'../ConfigurationFiles/02-Specific_ions.xlsx',
                   'score_cfg': r'../ConfigurationFiles/03-Score_weight_TG.xlsx'}

    usr_test_dct = {}
    usr_test_dct_keys = []
    mz_range = [600, 1000]  # default
    rt_range = [6, 10]  # default
    for usr_test in usr_test_lst:
        _test_dct = {'rank_score_filter': 27.5, 'score_filter': 27.5, 'isotope_score_filter': 75.0, 'ms_max': 0,
                     'pr_window': 0.75, 'ms2_infopeak_threshold': 0.001, 'ms2_hginfopeak_threshold': 0.001}
        if usr_test[0] in ['PC', 'PE', 'PA', 'PG', 'PI', 'PS', 'PIP']:
            lipid_class = usr_test[0]
            if lipid_class == 'PC':
                charge = '[M+HCOO]-'
            else:
                charge = '[M-H]-'
            vendor = usr_test[1]
            if vendor == 'waters':
                mzml = pl_mzml_waters
                mz_range = [650, 950]
                rt_range = [24, 27]  # max [24, 27]
            else:
                mzml = False
                pass

            _test_dct.update(pl_base_dct)

        elif usr_test[0] in ['TG']:
            lipid_class = usr_test[0]
            vendor = usr_test[1]
            charge = usr_test[2]
            if vendor == 'waters':
                mzml = tg_mzml_waters
                mz_range = [800, 1000]
                rt_range = [9, 15]  # max [9, 15]
            elif vendor == 'thermo':
                mzml = tg_mzml_thermo
                mz_range = [600, 1000]  # 850, 859
                rt_range = [20, 28]  # 25.6, 25.8
            elif vendor == 'sciex':
                mzml = tg_mzml_SCIEXS
                mz_range = [900, 1000]
                rt_range = [8, 13]
            elif vendor == 'agilent':
                mzml = tg_mzml_agilent
                mz_range = [600, 1000]
                rt_range = [9, 13]
            else:
                mzml = False

            _test_dct.update(tg_base_dct)
        elif usr_test[0] in ['DG']:
            # TODO (georgia.angelidou@uni-leipzig.de): remember to change the code below moer specific for DG
            lipid_class = usr_test[0]
            vendor = usr_test[1]
            charge = usr_test[2]

            if vendor == 'waters':
                mzml = tg_mzml_waters
                mz_range = [400, 800]
                rt_range = [6, 10]
            elif vendor == 'thermo':
                mzml = tg_mzml_thermo
                mz_range = [400, 800]
                rt_range = [6, 10]
            else:
                mzml = False

            _test_dct.update(tg_base_dct)
        else:

            lipid_class = False
            charge = False
            vendor = False
            mzml = False

        if mzml is not False and charge is not False:

            _cfg_dct = {'lipid_type': lipid_class, 'charge_mode': charge, 'vendor': vendor, 'mzml_path_str': mzml,
                        'rt_start': rt_range[0], 'rt_end': rt_range[1], 'mz_start': mz_range[0], 'mz_end': mz_range[1]}

            if vendor == 'waters':
                _cfg_dct['ms_ppm'] = 20
                _cfg_dct['ms2_ppm'] = 50
                _cfg_dct['hg_ppm'] = 200
                _cfg_dct['ms_th'] = 750
                _cfg_dct['ms2_th'] = 10
                _cfg_dct['hg_th'] = 10
                _cfg_dct['dda_top'] = 6

                _test_dct.update(_cfg_dct)
                usr_test_dct[usr_test[3]] = _test_dct
                usr_test_dct_keys.append(usr_test[3])

            elif vendor == 'thermo':
                _cfg_dct['ms_ppm'] = 5
                _cfg_dct['ms2_ppm'] = 20
                _cfg_dct['hg_ppm'] = 50
                _cfg_dct['ms_th'] = 10000
                _cfg_dct['ms2_th'] = 2000
                _cfg_dct['hg_th'] = 2000
                _cfg_dct['dda_top'] = 10

                _test_dct.update(_cfg_dct)
                usr_test_dct[usr_test[3]] = _test_dct
                usr_test_dct_keys.append(usr_test[3])
            elif vendor == 'sciex':
                _cfg_dct['ms_ppm'] = 10
                _cfg_dct['ms2_ppm'] = 60
                _cfg_dct['hg_ppm'] = 60
                _cfg_dct['ms_th'] = 5000
                _cfg_dct['ms2_th'] = 100  # Can get 2000/1000/750/500 depends how strict should be the identification
                _cfg_dct['hg_th'] = 1000
                _cfg_dct['dda_top'] = 15
                ms_ppm_SCIEXS = 10

                _test_dct.update(_cfg_dct)
                usr_test_dct[usr_test[3]] = _test_dct
                usr_test_dct_keys.append(usr_test[3])
            elif vendor == 'agilent':
                _cfg_dct['ms_ppm'] = 50
                _cfg_dct['ms2_ppm'] = 100
                _cfg_dct['hg_ppm'] = 100
                _cfg_dct['ms_th'] = 1000
                _cfg_dct['ms2_th'] = 10  # Can get 2000/1000/750/500 depends how strict should be the identification
                _cfg_dct['hg_th'] = 10
                _cfg_dct['dda_top'] = 4

                _test_dct.update(_cfg_dct)
                usr_test_dct[usr_test[3]] = _test_dct
                usr_test_dct_keys.append(usr_test[3])
        else:
            pass

    log_lst = []

    t0 = time.time()

    t_sum_lst = []

    # automatic identify the LipidHunter folder
    hunter_folder = os.path.dirname(os.getcwd())
    hunter_file_path = os.path.join(hunter_folder, 'LipidHunter.py')

    if os.path.isfile(hunter_file_path):
        print('\nLipidHunter folder', hunter_folder, '\n')
        print(usr_test_dct, '\n')

        for test_key in usr_test_dct_keys:
            if test_key in list(usr_test_dct.keys()):
                test_dct = usr_test_dct[test_key]
                t_str = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
                lipid_class = test_dct['lipid_type']

                cfg_dct = {'img_output_folder_str': r'../Test/results/%s_%s' % (test_key, t_str),
                           'xlsx_output_path_str': r'../Test/results/%s_%s.xlsx' % (test_key, t_str),
                           'hunter_folder': hunter_folder, 'img_type': u'png', 'img_dpi': 300,
                           'hunter_start_time': t_str, 'experiment_mode': 'LC-MS', 'rank_score': True,
                           'fast_isotope': False, 'core_number': core_count, 'max_ram': max_ram, 'tag_all_sn': True}

                test_dct.update(cfg_dct)

                print('>>>>>>>>>>>>>>>> START TEST: %s' % test_key)

                t, log_lst, export_df = huntlipids(test_dct, log_lst, save_fig=save_images)
                if t is not False:
                    print('>>>>>>>>>>>>>>>> TEST PASSED: %s in %.3f Sec <<<<<<<<<<<<<<<<\n' % (test_key, t))
                    t_sum_lst.append((test_key, 'PASSED', '%.3f Sec' % t, 'Identified: %i' % export_df.shape[0]))
                else:
                    print('>>>>>>>>!!!!!!!! TEST FAILED: %s !!!!!!!<<<<<<<<\n' % test_key)
                    t_sum_lst.append((test_key, 'FAILED', '', 'Identified: 0'))

    else:
        print('!!! Invalid LipidHunter folder', hunter_folder)
        for test_key in usr_test_dct_keys:
            print('>>>>>>>>!!!!!!!! TEST FAILED: %s !!!!!!!<<<<<<<<\n' % test_key)
            t_sum_lst.append((test_key, 'FAILED', '', 'Identified: 0'))

    t_end = time.time() - t0

    print('Test run in plan: ', ', '.join(usr_test_dct_keys))
    print('With Max Core = %i and RAM = %i GB' % (core_count, max_ram))
    if len(t_sum_lst) > 0:
        for t_info in t_sum_lst:
            print('    '.join(t_info))
    print('\n=============== ALL TEST FINISHED in %.3f Sec ===============' % t_end)
