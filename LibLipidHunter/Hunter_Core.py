# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
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
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
#     Developer Georgia Angelidou georgia.angelidou@uni-leipzig.de

import math
import multiprocessing
from multiprocessing import Pool
import os
import sys
from sys import platform
import time
from typing import Tuple, Union

from numpy import int64
import pandas as pd

from LibLipidHunter.LipidComposer import LipidComposer
from LibLipidHunter.SpectraReader import extract_mzml
from LibLipidHunter.SpectraReader import get_spectra
from LibLipidHunter.SpectraReader import get_xic_from_pl
from LibLipidHunter.SpectraReader import get_spec_info
from LibLipidHunter.LogPageCreator import LogPageCreator
from LibLipidHunter.PrecursorHunter import PrecursorHunter
from LibLipidHunter.ScoreHunter import get_lipid_info
from LibLipidHunter.PanelPlotter import gen_plot
from LibLipidHunter.HuntManager import save_hunt
from LibLipidHunter.HuntManager import gen_html_report


def save_output(output_path: str, output_df: pd.DataFrame, output_name: str = 'output'):
    """
    function to save output in a safer way with few checks.

    Args:
        output_path (str): The path to save the output file.
        output_df (pd.DataFrame): The DataFrame to be saved.
        output_name (str, optional): The name of the DataFrame to be displayed in debug information.

    """

    output_directory = os.path.dirname(output_path)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    if output_path.endswith('.xlsx'):
        try:
            output_df.to_excel(output_path, index=False)
            print(f'[OUTPUT] ==> Prepare to save {output_name} as: ', output_path)
        except IOError:
            output_df.to_excel('{0}-{1}{2}'.format(output_path[:-5], int(time.time()), '.xlsx'), index=False)
            print(output_path)
    elif output_path.endswith('.csv'):
        try:
            output_df.to_csv(output_path, index=False)
            print(f'[OUTPUT] ==> Prepare to save {output_name} as: ', output_path)
        except IOError:
            output_df.to_csv('{0}-{1}{2}'.format(output_path[:-5], int(time.time()), '.csv'), index=False)
            print(output_path)
    else:
        try:
            output_df.to_excel(f'{output_path}.xlsx', index=False)
            print(f'[OUTPUT] ==> Prepare to save {output_name} as: ', f'{output_path}.xlsx')
        except IOError:
            output_df.to_excel('{0}-{1}{2}'.format(output_path[:-5], int(time.time()), '.xlsx'), index=False)
            print(output_path)
    print('[OUTPUT] ==> {output_name} saved ...')


def huntlipids(param_dct: dict, error_lst: list,
               save_fig: bool = True, save_session: bool = False) \
        -> Tuple[Union[float, bool], list, Union[pd.DataFrame, bool]]:
    """
    This is the core workflow of lipidhunter.
    The cmd_lipidhunter and LipidHunter GUI transfer the parsed Lipidhunter configurations as param_dct.
    This function read all parameters and control the whole workflow.

    Examples:
        {'fawhitelist_path_str': r'D:\lipidhunter\ConfigurationFiles\FA_Whitelist.xlsx',
        'mzml_path_str': r'D:\lipidhunter\test\mzML\PL_Neg_Waters_qTOF.mzML',
        'img_output_folder_str': r'D:\lipidhunter\Temp\Test2',
        'xlsx_output_path_str': r'D:\lipidhunter\Temp\Test2\t2.xlsx',
        'lipid_specific_cfg': r'D:\lipidhunter\ConfigurationFiles\PL_specific_ion_cfg.xlsx',
        'hunter_start_time': '2017-12-21_15-27-49',
        'vendor': 'waters', 'experiment_mode': 'LC-MS', 'lipid_class': 'PC', 'charge_mode': '[M+HCOO]-',
        'rt_start': 20.0, 'rt_end': 25.0, 'mz_start': 700.0, 'mz_end': 800.0,
        'rank_score': True, 'rank_score_filter': 27.5, 'score_filter': 27.5,
        'isotope_score_filter': 75.0, 'fast_isotope': False,
        'ms_th': 1000, 'ms_ppm': 20, 'ms_max': 0, 'pr_window': 0.75, 'dda_top': 6,
        'ms2_th': 10, 'ms2_ppm': 50, 'ms2_infopeak_threshold': 0.001,
        'score_cfg': r'D:\lipidhunter\ConfigurationFiles\Score_cfg.xlsx',
        'hunter_folder': r'D:\lipidhunter',
        'core_number': 3, 'max_ram': 5, 'img_type': u'png', 'img_dpi': 300}

    Args:
        param_dct (dict): The overall parameters received for processing.
        error_lst (list): All errors collected before this function.
        save_fig (bool): Can be set to False to skip image generation (not recommended).
        save_session (bool): Can be set to True to save session as a file using pickle.

    Returns:
        tot_run_time (float): The total run time as float or False if error occurs.
        error_lst (list): List of errors collected during the run.
        output_df (pd.DataFrame, bool): The output DataFame or False if error occurs.

    """

    print('[INFO] --> Hunter Core Start...')

    # initialization of all parameters

    start_time = time.clock()
    lipidcomposer = LipidComposer()

    usr_lipid_class = param_dct['lipid_class']
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
    usr_ms1_ppm = param_dct['ms_ppm']
    usr_ms2_ppm = param_dct['ms2_ppm']
    usr_ms1_precision = usr_ms1_ppm * 1e-6
    usr_ms2_precision = usr_ms2_ppm * 1e-6
    # usr_rank_score_filter = param_dct['rank_score_filter']
    # usr_score_filter = param_dct['score_filter']
    # usr_isotope_score_filter = param_dct['isotope_score_filter']
    # usr_ms2_info_th = param_dct['ms2_infopeak_threshold']
    # usr_rank_mode = param_dct['rank_score']
    # usr_fast_isotope = param_dct['fast_isotope']
    if 'xic_ppm' in list(param_dct.keys()):
        usr_xic_ppm = int(param_dct['xic_ppm'])
    else:
        usr_xic_ppm = usr_ms1_ppm

    # parameters from settings tab
    usr_core_num = param_dct['core_number']
    usr_max_ram = param_dct['max_ram']

    # usr_dpi = param_dct['img_dpi']
    # usr_img_type = param_dct['img_type']

    hunter_start_time_str = param_dct['hunter_start_time']

    try:
        if param_dct['debug_mode'] == 'ON':
            save_session = True
        else:
            pass
    except (KeyError, AttributeError):
        print('Debug mode == off')

    # Check platform for multiprocessing settings

    if platform == "linux" or platform == "linux2":  # linux
        if usr_core_num > 1:
            os_typ = 'linux_multi'
            print('[INFO] --> LipidHunter Running on >>> Linux with multiprocessing mode ...')
        else:
            os_typ = 'linux_single'
            print('[INFO] --> LipidHunter Running on >>> Linux with single core mode ...')
    elif platform == "win32":  # Windows
        os_typ = 'windows'
        print('[INFO] --> LipidHunter Running on >>> Windows ...')
    elif platform == "darwin":  # macOS
        if usr_core_num > 1:
            os_typ = 'linux_multi'
            print('[INFO] --> LipidHunter Running on >>> macOS with multiprocessing mode ...')
        else:
            os_typ = 'linux_single'
            print('[INFO] --> LipidHunter Running on >>> macOS with single core mode ...')
    else:
        usr_core_num = 1
        param_dct['core_number'] = 1
        os_typ = 'linux_single'
        print('[INFO] --> LipidHunter Running on >>> unoptimized system ... %s' % platform)
        print('[WARNING] !!! Force to use single core mode !!')

    print('[INFO] --> Start to process >>>')
    print('[INFO] --> Lipid class: %s >>>' % usr_lipid_class)

    composer_param_dct = {'fa_whitelist': usr_fa_xlsx, 'lipid_class': usr_lipid_class,
                          'charge_mode': usr_charge, 'exact_position': 'FALSE'}

    existed_lipid_master_path = ''  # reserved to load existing lipid master table for later development
    use_existed_lipid_master = False  # reserved to re-use lipid master table
    save_lipid_master_table = False
    if 'debug_mode' in list(param_dct.keys()):
        if param_dct['debug_mode'] == 'ON':
            if 'lipid_master_table' in list(param_dct.keys()):
                existed_lipid_master_path = param_dct['lipid_master_table']
                if os.path.isfile(existed_lipid_master_path):
                    use_existed_lipid_master = True
                else:
                    print('[ERROR] !!! Failed to load existed Lipid Master table: %s', existed_lipid_master_path)
            else:
                save_lipid_master_table = True

        # reserved to re-use lipid master table
        # if 'save_lipid_master_table' in list(param_dct.keys()):
        #     if param_dct['save_lipid_master_table'] == 'CSV':
        #         save_lipid_master_table = True

    if use_existed_lipid_master is False:
        try:
            print('[INFO] --> Start to generate Lipid Master Table ...')
            t_lm_0 = time.time()
            usr_lipid_master_df = lipidcomposer.compose_lipid(param_dct=composer_param_dct, ms2_ppm=usr_ms2_ppm)
            print('[INFO] --> Lipid Master Table generated >>> in %.2f sec' % (time.time() - t_lm_0))

        except Exception as e:
            print('[ERROR] !!! Failed to predict Lipid structures...', e)
            error_lst.append(e)
            error_lst.append('[ERROR] !!! Some files missing...')
            error_lst.append('... ... Please check your settings in the configuration file ...')
            return False, error_lst, False
    else:
        try:
            print('[INFO] --> Try to use existed Lipid Master table: %s' % existed_lipid_master_path)
            usr_lipid_master_df = pd.read_csv(existed_lipid_master_path)
            print('[INFO] --> Lipid Master table loaded >>>', usr_lipid_master_df.shape[0])
        except Exception as e:
            print(e)
            error_lst.append(e)
            return False, error_lst, False

    if isinstance(usr_lipid_master_df, pd.DataFrame):

        if not usr_lipid_master_df.empty:
            pass
        else:
            print('[ERROR] !!! Failed to generate LipidMaster Table...')
            error_lst.append('[ERROR] !!! Failed to generate LipidMaster Table...')
            error_lst.append('... ... Please check if Lipid Class and FA are marked in FA whitelist...')
            return False, error_lst, False
    else:
        print('[ERROR] !!! Failed to generate LipidMaster Table...')
        error_lst.append('[ERROR] !!! Failed to generate LipidMaster Table...')
        error_lst.append('... ... Please check if Lipid Class and FA are marked in FA whitelist...')
        return False, error_lst, False

    if save_lipid_master_table is True:
        log_master_name = 'Lipid_Master_%s.csv' % hunter_start_time_str
        log_master_name = os.path.join(output_folder, log_master_name)
        print('[OUTPUT] ==> Lipid Master table Saved as: ', log_master_name)
        save_output(log_master_name, usr_lipid_master_df, output_name='Lipid Master table')
    else:
        pass

    # for TG has the fragment of neutral loss of the FA and the fragments for the MG
    usr_fa_df = lipidcomposer.calc_fa_query(usr_lipid_class, usr_fa_xlsx, ms2_ppm=usr_ms2_ppm)
    # del lipidcomposer
    if usr_fa_df is False:
        print('[ERROR] !!! Failed to generate FA info table ...\n')
        error_lst.append('[ERROR] !!! Failed to generate FA info table ...\n')
        return False, error_lst, False

    lipid_info_df = usr_lipid_master_df  # make a copy to avoid pandas warning when modify it later

    # cut lib info to the user defined m/z range
    # TODO (georgia.angelidou@uni-leipzig.de): support for the sphingomyelins and ceramides
    pos_charge_lst = ['[M+H]+', '[M+Na]+', '[M+NH4]+']
    neg_charge_lst = ['[M-H]-', '[M+HCOO]-', '[M+CH3COO]-']
    if usr_charge in neg_charge_lst:
        if usr_lipid_class in ['PC', 'LPC']:
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
        if usr_lipid_class == 'DG':
            if usr_charge == '[M+NH4]+':
                lipid_info_df = lipid_info_df[
                    (mz_start <= lipid_info_df['[M+NH4]+_MZ']) & (lipid_info_df['[M+NH4]+_MZ'] <= mz_end)]
            elif usr_charge == '[M+H]+':
                lipid_info_df = lipid_info_df[
                    (mz_start <= lipid_info_df['[M+H]+_MZ']) & (lipid_info_df['[M+H]+_MZ'] <= mz_end)]
    else:
        error_lst.append('Lipid class or charge NOT supported.  User input lipid class = %s, charge = %s. '
                         % (usr_lipid_class, usr_charge))

    # TODO(zhixu.ni@uni-leipzig.de): Add more error to the error_lst.

    pr_hunter = PrecursorHunter(lipid_info_df, param_dct, os_type=os_typ)

    output_df = pd.DataFrame()

    print('[INFO] --> Start to process data ...')
    print('[INFO] --> Lipid class: %s' % usr_lipid_class)

    # generate the Weight factor df
    usr_weight_df = pd.read_excel(score_cfg, index_col='Type')

    print('[INFO] --> Start to parse mzML')

    # Preparation finished.
    # Lipid identification workflow start.

    # Extract all spectra from mzML to pandas DataFrame
    usr_scan_info_df, usr_spectra_pl, ms1_xic_df = extract_mzml(usr_mzml, usr_rt_range, dda_top=usr_dda_top,
                                                                ms1_threshold=usr_ms1_threshold,
                                                                ms2_threshold=usr_ms2_threshold,
                                                                ms1_precision=usr_ms1_precision,
                                                                ms2_precision=usr_ms2_precision,
                                                                vendor=usr_vendor, ms1_max=usr_ms1_max)

    print('[INFO] --> MS1_XIC_df.shape', ms1_xic_df.shape)
    # Find all possible precursor according to lipid master table
    ms1_obs_pr_df = pr_hunter.get_matched_pr(usr_scan_info_df, usr_spectra_pl, ms1_max=usr_ms1_max,
                                             core_num=usr_core_num, max_ram=usr_max_ram)

    if ms1_obs_pr_df is False:
        print('[WARNING] !!! NO suitable precursor --> Check settings!!\n')
        error_lst.append('[WARNING] !! NO suitable precursor --> Check settings!!\n')
        return False, error_lst, False

    print('[INFO] --> ms1 precursor matched')

    # Remove bad precursors, keep the matched scans by DDA_rank and scan number
    # Build unique str identifier for each scan with scan_number00dda_rank use numpy.int64 to avoid large scan_number
    # The unique identifier is important to index scans
    usr_scan_info_df['scan_checker'] = (usr_scan_info_df['scan_number'].astype(int64).astype(str).str.
                                        cat(usr_scan_info_df['DDA_rank'].astype(int64).astype(str), sep='_'))
    ms1_obs_pr_df['scan_checker'] = (ms1_obs_pr_df['scan_number'].astype(int64).astype(str).str.
                                     cat(ms1_obs_pr_df['DDA_rank'].astype(int64).astype(str).astype(str), sep='_'))

    checked_info_df = ms1_obs_pr_df[ms1_obs_pr_df['scan_checker'].isin(usr_scan_info_df['scan_checker']
                                                                       .tolist())].copy()
    checked_info_df.sort_values(by=['scan_checker', 'Lib_mz'], ascending=[True, True], inplace=True)

    # # reserved for further debug_mode options
    # if 'debug_mode' in list(param_dct.keys()):
    #     if param_dct['debug_mode'] == 'ON':
    #         usr_scan_info_df.to_csv(os.path.join(output_folder, 'usr_scan_info.csv'))
    #         ms1_obs_pr_df.to_csv(os.path.join(output_folder, 'ms1_obs_pr_df.csv'))
    #         checked_info_df.to_csv(os.path.join(output_folder, 'checked_info_df.csv'))

    if checked_info_df.empty:
        print('[ERROR] !!! No identification in pre-match steps !!')
        error_lst.append('!! No identification in pre-match steps !!\n')
        return False, error_lst, False
    else:
        print('[INFO] --> features identified in the pre-match: ', checked_info_df.shape[0])

    ms1_xic_mz_lst = sorted(set(ms1_obs_pr_df['MS1_XIC_mz'].values.tolist()))
    print('ms1_xic_mz_lst', len(ms1_xic_mz_lst))
    print(ms1_xic_mz_lst)
    print('[INFO] --> Start to extract XIC')

    # Distribute task according to core numbers
    if len(ms1_xic_mz_lst) >= 3 * usr_core_num:
        sub_len = int(math.ceil(len(ms1_xic_mz_lst) / usr_core_num))
        core_key_list = [ms1_xic_mz_lst[k: k + sub_len] for k in range(0, len(ms1_xic_mz_lst), sub_len)]
    else:
        core_key_list = [ms1_xic_mz_lst]

    # Start multiprocessing to get XIC according to precursor m/z found above
    print('[STATUS] >>> Start multiprocessing to get XIC ==> ==> ==> Number of Cores: %i' % usr_core_num)
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
                    print('[STATUS] >>> Core #%i ==> ...... processing ......' % worker_count)
                    print(core_list)
                    # TODO (georgia.angelidou@uni-leipzig.de): maybe can be combine or stay like this????
                    xic_result = parallel_pool.apply_async(get_xic_from_pl, args=(core_list, ms1_xic_df, usr_xic_ppm,
                                                                                  os_typ, queue))
                    worker_count += 1
                    xic_results_lst.append(xic_result)

            parallel_pool.close()
            parallel_pool.join()
            # del ms1_xic_df
            # TODO (georgia.angelidou@uni-leipzig.de): can be done before the join??????
            for xic_result in xic_results_lst:
                try:
                    sub_xic_dct = xic_result.get()
                    if len(list(sub_xic_dct.keys())) > 0:
                        xic_dct.update(sub_xic_dct)

                except (KeyError, SystemError, ValueError):
                    pass
            del xic_results_lst

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
                    print('[STATUS] >>> Core #%i ==> ...... processing ......' % worker_count)
                    print(core_list)
                    job = multiprocessing.Process(target=get_xic_from_pl, args=(core_list, ms1_xic_df, usr_xic_ppm,
                                                                                os_typ, queue))
                    worker_count += 1
                    jobs.append(job)
                    job.start()
                    xic_results_lst.append(queue.get())
            del ms1_xic_df

            for j in jobs:
                j.join()

            # TODO (georgia,angelidou@uni-leipzig.de):
            #  question does these part can be move before the different jobs are join?????
            for xic_result in xic_results_lst:
                try:
                    if len(list(xic_result.keys())) > 0:
                        xic_dct.update(xic_result)
                except (KeyError, SystemError, ValueError):
                    pass
            del xic_results_lst
    else:
        print('[INFO] --> Using single core mode...')
        queue = ''
        worker_count = 1
        for core_list in core_key_list:
            if isinstance(core_list, tuple) or isinstance(core_list, list):
                if None in core_list:
                    core_list = [x for x in core_list if x is not None]
                else:
                    pass
                print('[STATUS] >>> Core #1 Part %i ==> ...... processing ......' % worker_count)
                print(core_list)
                sub_xic_dct = get_xic_from_pl(core_list, ms1_xic_df, usr_xic_ppm, os_typ, queue)
                worker_count += 1

                if len(list(sub_xic_dct.keys())) > 0:
                    xic_dct.update(sub_xic_dct)
        del ms1_xic_df

    if len(list(xic_dct.keys())) == 0:
        print('[ERROR] !!! No precursor for XIC found !!')
        error_lst.append('!! No precursor for XIC found !!\n')
        return False, error_lst, False
    else:
        print('[INFO] --> Number of XIC extracted: %i' % len(list(xic_dct.keys())))

    print('[STATUS] >>> Start to Hunt for Lipids !!')
    # All XIC extracted
    # Start to process further identification and scoring
    # Stat to use get_spec_info to match MS2 to precursor and XIC with multiprocessing

    # Todo(zhixu.ni@uni-leipzig.de): extract the multiprocessing part into function

    target_ident_lst = []
    lipid_spec_info_dct = {}
    spec_part_key_lst = []
    split_seg = 1
    lipid_all_scan_num = checked_info_df['scan_checker'].nunique()

    # Checked_info_df may become very large for TG in Thermo files, need to be divided for Multiprocessing
    checked_info_df_groups = checked_info_df.groupby(['Formula', 'scan_checker'])
    checked_info_key_lst = list(checked_info_df_groups.groups.keys())
    checked_info_key_num = len(checked_info_key_lst)
    print('[INFO] --> Total number of spectra features: %i | Number of proposed Formula: %i'
          % (lipid_all_scan_num, checked_info_key_num))

    if lipid_all_scan_num > (usr_core_num * 100):

        # Split tasks into few parts to avoid core waiting in multiprocessing
        if usr_core_num * 200 < lipid_all_scan_num <= usr_core_num * 400:
            split_seg = 2
        elif usr_core_num * 400 < lipid_all_scan_num:
            split_seg = 3
        else:
            split_seg = 2

        spec_part_len = int(math.ceil(lipid_all_scan_num / split_seg))  # Number of scan in each part
        if spec_part_len > 200:
            spec_part_len = 200

        spec_part_lst = [list(checked_info_df['scan_checker'].unique())[k: k + spec_part_len]
                         for k in range(0, lipid_all_scan_num, spec_part_len)]

        for part_lst in spec_part_lst:
            if None in part_lst:
                part_lst = [x for x in part_lst if x is not None]
            spec_sub_len = int(math.ceil(len(part_lst) / usr_core_num))
            scan_sub_lst = [part_lst[k: k + spec_sub_len] for k in range(0, len(part_lst), spec_sub_len)]
            for scan_lst in scan_sub_lst:
                _tmp_lipid_info_df = checked_info_df[checked_info_df['scan_checker'].isin(scan_lst)]
                _tmp_info_groups = _tmp_lipid_info_df.groupby(['Formula', 'scan_checker'])
                _tmp_group_key_lst = list(_tmp_info_groups.groups.keys())
                _tmp_group_key_num = len(_tmp_group_key_lst)
                # chk_info_df_lst.append([_tmp_lipid_info_df, _tmp_info_groups])
                spec_sub_key_lst = [_tmp_group_key_lst[k: k + spec_sub_len]
                                    for k in range(0, _tmp_group_key_num, spec_sub_len)]
                spec_part_key_lst.append((spec_sub_key_lst, _tmp_info_groups))

        print(
            '[INFO] --> Distributed tasks into %i parts with max %i features ...' % (len(spec_part_lst), spec_part_len))

    else:

        spec_sub_len = int(math.ceil(checked_info_key_num / usr_core_num))
        spec_sub_key_lst = [checked_info_key_lst[k: k + spec_sub_len] for k in
                            range(0, checked_info_key_num, spec_sub_len)]
        spec_part_key_lst.append((spec_sub_key_lst, checked_info_df_groups))
    del checked_info_df_groups

    if usr_core_num > 1:
        part_tot = len(spec_part_key_lst)
        part_counter = 1
        queue = ''

        for spec_sub_lst in spec_part_key_lst:

            spec_sub_key_lst = spec_sub_lst[0]
            sub_info_groups = spec_sub_lst[1]

            if part_tot == 1:
                print('[STATUS] >>> Start multiprocessing to get Spectra info ==> Max Number of Cores: %i'
                      % usr_core_num)
            else:
                print('[STATUS] >>> Start multiprocessing to get Spectra info ==> Part %i / %i '
                      '--> Max Number of Cores: %i'
                      % (part_counter, part_tot, usr_core_num))

            spec_results_lst = []

            if os_typ == 'windows':
                parallel_pool = Pool(usr_core_num)
                # usr_queue = ''
                worker_count = 1
                for _sub_lst in spec_sub_key_lst:
                    if isinstance(_sub_lst, tuple) or isinstance(_sub_lst, list):
                        if None in _sub_lst:
                            _sub_lst = [x for x in _sub_lst if x is not None]
                        else:
                            pass
                        print('[STATUS] >>> Core #%i ==> ...... processing ......' % worker_count)
                        spec_result = parallel_pool.apply_async(get_spec_info, args=(_sub_lst, sub_info_groups,
                                                                                     usr_scan_info_df, os_typ, queue))
                        worker_count += 1
                        if worker_count > usr_core_num:
                            worker_count = 1
                        spec_results_lst.append(spec_result)

                parallel_pool.close()
                parallel_pool.join()

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
                        print('[STATUS] >>> Core #%i ==> ...... processing ......' % worker_count)
                        job = multiprocessing.Process(target=get_spec_info, args=(_sub_lst, sub_info_groups,
                                                                                  usr_scan_info_df, os_typ, queue_spec))
                        worker_count += 1
                        jobs.append(job)
                        job.start()
                        spec_results_lst.append(queue_spec.get())
                for j in jobs:
                    j.join()
            del sub_info_groups
            #  Merge multiprocessing results
            for spec_result in spec_results_lst:

                if os_typ == 'windows':
                    try:
                        sub_spec_dct = spec_result.get()
                        # print(sub_spec_dct)
                        if len(list(sub_spec_dct.keys())) > 0:
                            lipid_spec_info_dct.update(sub_spec_dct)
                    except (KeyError, SystemError, ValueError):
                        print('[ValueError] !!! must supply a tuple to get_group with multiple grouping keys ...')
                else:  # for linux
                    try:
                        if len(list(spec_result.keys())) > 0:
                            lipid_spec_info_dct.update(spec_result)
                    except (KeyError, SystemError, ValueError):
                        print('[ValueError] !!! must supply a tuple to get_group with multiple grouping keys ...')
            # del spec_result
            # del spec_results_lst
            if part_tot == 1:
                print('[STATUS] >>> multiprocessing results merged')
            else:
                print('[STATUS] >>> multiprocessing results merged ==> Part %i / %i '
                      % (part_counter, part_tot))
            part_counter += 1
    else:
        print('[INFO] --> Using single core mode...')
        queue = ''
        worker_count = 1
        for spec_sub_lst in spec_part_key_lst:
            spec_sub_key_lst = spec_sub_lst[0]
            sub_info_groups = spec_sub_lst[1]
            for _sub_lst in spec_sub_key_lst:
                if isinstance(_sub_lst, tuple) or isinstance(_sub_lst, list):
                    if None in _sub_lst:
                        _sub_lst = [x for x in _sub_lst if x is not None]
                    else:
                        if isinstance(_sub_lst[0], float):
                            _sub_lst3 = ()
                            _sub_lst3 = _sub_lst3 + (_sub_lst,)
                            _sub_lst = _sub_lst3

                    print('[STATUS] >>> Core #%i ==> ...... processing ......' % worker_count)
                    sub_spec_dct = get_spec_info(_sub_lst, sub_info_groups, usr_scan_info_df, os_typ, queue)

                    if len(list(sub_spec_dct.keys())) > 0:
                        lipid_spec_info_dct.update(sub_spec_dct)

                else:
                    pass

    print('lipid_spec_info_dct', len(list(lipid_spec_info_dct.keys())))

    # Single process ONLY. usr_spectra_pl is too big in RAM --> RAM leaking during copy
    lipid_spec_dct = {}
    # TODO (georgia.angelidou@uni-leipzig.de): when for loop change below line will be remove
    spec_info_key_lst = list(lipid_spec_info_dct.keys())

    # TODO (georgia.angelidou@uni-leipzig.de):
    #  Need to change since it is used only once is not necessary to introduce a new parameter
    #  change to :
    #  for _spec_group_key in list(lipid_spec_info_dct.keys())

    for _spec_group_key in spec_info_key_lst:
        _spec_info_dct = lipid_spec_info_dct[_spec_group_key]
        _usr_ms2_pr_mz = _spec_info_dct['MS2_PR_mz']
        _usr_ms2_dda_rank = _spec_info_dct['DDA_rank']
        _usr_ms2_scan_id = _spec_info_dct['scan_number']
        _usr_mz_lib = _spec_info_dct['Lib_mz']

        usr_spec_info_dct = get_spectra(_usr_ms2_pr_mz, _usr_mz_lib, _usr_ms2_dda_rank, _usr_ms2_scan_id,
                                        ms1_xic_mz_lst, usr_scan_info_df, usr_spectra_pl,
                                        dda_top=usr_dda_top, ms1_precision=usr_ms1_precision, vendor=usr_vendor)
        lipid_spec_dct[_spec_group_key] = usr_spec_info_dct

    found_spec_key_lst = list(lipid_spec_dct.keys())
    found_spec_key_lst = sorted(found_spec_key_lst, key=lambda x: x[0])
    spec_key_num = len(found_spec_key_lst)

    # Peak info matched
    # Parse specific peak info
    # Lipid class specific triggers

    pl_class_lst = ['PA', 'PC', 'PE', 'PG', 'PI', 'PS', 'PIP']
    lpl_class_lst = ['LPA', 'LPC', 'LPE', 'LPG', 'LPI', 'LPS', 'LPIP']
    pl_neg_chg_lst = ['[M-H]-', '[M+HCOO]-', '[M+CH3COO]-']
    tg_class_lst = ['TG', 'DG']
    tg_pos_chg_lst = ['[M+NH4]+', '[M+H]+', '[M+Na]+']
    # Todo(zhixu.ni@uni-leipzig.de): Other classes such as SM, Cer and CL here.

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
    elif usr_lipid_class in lpl_class_lst and usr_charge in pl_neg_chg_lst:
        charge_mode = 'NEG'
        usr_key_frag_df = pd.read_excel(key_frag_cfg)
        usr_key_frag_df = usr_key_frag_df.query('EXACTMASS > 0')
        # get the information from the following columns
        usr_key_frag_df = usr_key_frag_df[['CLASS', 'TYPE', 'EXACTMASS', 'PR_CHARGE', 'LABEL', 'CHARGE_MODE']]
        # find key peaks for the target PL class
        target_frag_df = usr_key_frag_df.query(r'CLASS == "%s" and TYPE == "FRAG" and PR_CHARGE == "%s"'
                                               % (usr_lipid_class[1:], usr_charge))
        target_nl_df = usr_key_frag_df.query(r'CLASS == "%s" and TYPE == "NL" and PR_CHARGE == "%s"'
                                             % (usr_lipid_class[1:], usr_charge))
        # add precursor to the list
        target_pr_df = pd.DataFrame(data={'CLASS': usr_lipid_class, 'TYPE': 'NL', 'EXACTMASS': 0.0,
                                          'PR_CHARGE': usr_charge, 'LABEL': 'PR', 'CHARGE_MODE': 'NEG'}, index=['PR'])
        target_nl_df = target_nl_df.append(target_pr_df)
        target_nl_df.reset_index(drop=True, inplace=True)

        # extract info for other classes
        other_frag_df = usr_key_frag_df.query('CLASS != "%s" and TYPE == "FRAG" and CHARGE_MODE == "%s"'
                                              % (usr_lipid_class[1:], charge_mode))
        other_nl_df = usr_key_frag_df.query('CLASS != "%s" and TYPE == "NL" and CHARGE_MODE == "%s"'
                                            % (usr_lipid_class[1:], charge_mode))
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

    print('[INFO] --> Key FRAG Dict Generated ...')

    # Start to get rank score using get_lipid_info with multiprocessing

    lipid_info_img_lst = []
    lipid_part_key_lst = []
    split_seg = 1

    if spec_key_num > (usr_core_num * 24):

        # Split tasks into few parts to avoid core waiting in multiprocessing
        # Problem with the 1 Core run. Temporary solution
        if usr_core_num * 24 < spec_key_num <= usr_core_num * 48 and usr_core_num != 1:
            split_seg = 2
        elif usr_core_num * 48 < spec_key_num <= usr_core_num * 96 and usr_core_num != 1:
            split_seg = 3
        elif usr_core_num * 96 < spec_key_num and usr_core_num != 1:
            split_seg = 4
        else:
            split_seg = 1
        lipid_part_len = int(math.ceil(spec_key_num / split_seg))
        lipid_part_lst = [found_spec_key_lst[k: k + lipid_part_len] for k in range(0, spec_key_num,
                                                                                   lipid_part_len)]
        print('[INFO] --> lipid_part_number: ', len(lipid_part_lst), ' lipid_part_len:', lipid_part_len)

        for part_lst in lipid_part_lst:
            if None in part_lst:
                part_lst = [x for x in part_lst if x is not None]
            lipid_sub_len = int(math.ceil(len(part_lst) / usr_core_num))
            # print('lipid_sub_len', lipid_sub_len)
            pre_lipid_sub_key_lst = [part_lst[k: k + lipid_sub_len] for k in range(0, len(part_lst), lipid_sub_len)]
            lipid_sub_key_lst = []
            for _core_key_lst in pre_lipid_sub_key_lst:
                _core_key_df = pd.DataFrame(_core_key_lst, columns=['Formula', 'scan_checker'])
                _core_scan_lst = _core_key_df['scan_checker'].tolist()
                _core_chk_info_df = checked_info_df[checked_info_df['scan_checker'].isin(_core_scan_lst)]
                _core_chk_groups = _core_chk_info_df.groupby(['Formula', 'scan_checker'])

                lipid_sub_key_lst.append((_core_key_lst, _core_chk_info_df, _core_chk_groups))

            lipid_part_key_lst.append(lipid_sub_key_lst)

    else:
        lipid_sub_len = int(math.ceil(spec_key_num / usr_core_num))
        pre_lipid_sub_key_lst = [found_spec_key_lst[k: k + lipid_sub_len] for k in
                                 range(0, spec_key_num, lipid_sub_len)]
        lipid_sub_key_lst = []
        for _core_key_lst in pre_lipid_sub_key_lst:
            _core_key_df = pd.DataFrame(_core_key_lst, columns=['Formula', 'scan_checker'])
            _core_scan_lst = _core_key_df['scan_checker'].tolist()
            _core_chk_info_df = checked_info_df[checked_info_df['scan_checker'].isin(_core_scan_lst)]
            _core_chk_groups = _core_chk_info_df.groupby(['Formula', 'scan_checker'])

            # TODO (georgia.angelidou@uni-leipzig.de):  unnecesary decleration of new parameter
            #  It can directly be save in the destination parameter:
            #  lipid_part_key_lst.append((_core_key_lst, _core_chk_info_df, _core_chk_groups))
            lipid_sub_key_lst.append((_core_key_lst, _core_chk_info_df, _core_chk_groups))

        lipid_part_key_lst.append(lipid_sub_key_lst)

    part_tot = len(lipid_part_key_lst)
    part_counter = 1
    queue = ''

    if usr_core_num > 1:

        for lipid_sub_key_lst in lipid_part_key_lst:

            lipid_info_results_lst = []

            if part_tot == 1:
                print('[STATUS] >>> Start multiprocessing to get Score ==> Max Number of Cores: %i' % usr_core_num)
                try:
                    worker_feature_count = len(lipid_sub_key_lst[0][0])
                    print('[STATUS] >>> Start multiprocessing to get Score ==> Max Number of Cores: %i | '
                          'x%i Features each'
                          % (usr_core_num, worker_feature_count))
                except Exception as _e:
                    print('[STATUS] >>> Start multiprocessing to get Score ==> Max Number of Cores: %i' % usr_core_num)
                    print('[Exception] Can not get the number of features distributed to each core...', _e)
            else:
                try:
                    worker_feature_count = len(lipid_sub_key_lst[0][0])
                    print('[STATUS] >>> Start multiprocessing to get Score ==> Part %i / %i '
                          '--> Max Number of Cores: %i | x%i Features each'
                          % (part_counter, part_tot, usr_core_num, worker_feature_count))
                except Exception as _e:
                    print('[STATUS] >>> Start multiprocessing to get Score ==> Part %i / %i --> Max Number of Cores: %i'
                          % (part_counter, part_tot, usr_core_num))
                    print('[Exception] ... Can not get the number of features distributed to each core...', _e)

            if os_typ == 'windows':
                parallel_pool = Pool(usr_core_num)

                worker_count = 1
                for lipid_sub in lipid_sub_key_lst:
                    if isinstance(lipid_sub, tuple) or isinstance(lipid_sub, list):

                        lipid_sub_lst = lipid_sub[0]
                        _chk_info_df = lipid_sub[1]
                        _chk_info_gp = lipid_sub[2]

                        if None in lipid_sub_lst:
                            lipid_sub_lst = [x for x in lipid_sub_lst if x is not None]
                        else:
                            pass
                        if isinstance(lipid_sub_lst[0], tuple) or isinstance(lipid_sub_lst[0], list):
                            lipid_sub_dct = {k: lipid_spec_dct[k] for k in lipid_sub_lst}
                        else:
                            lipid_sub_dct = {lipid_sub_lst: lipid_spec_dct[lipid_sub_lst]}
                            lipid_sub_lst = tuple([lipid_sub_lst])
                        print('[STATUS] >>> Core #%i ==> ...... processing ......' % worker_count)
                        if len(list(lipid_sub_dct.keys())) > 0:
                            lipid_info_result = parallel_pool.apply_async(get_lipid_info,
                                                                          args=(param_dct, usr_fa_df, _chk_info_df,
                                                                                _chk_info_gp, lipid_sub_lst,
                                                                                usr_weight_df, key_frag_dct,
                                                                                lipid_sub_dct, xic_dct,
                                                                                worker_count,
                                                                                save_fig, os_typ, queue))

                            lipid_info_results_lst.append(lipid_info_result)
                            worker_count += 1

                parallel_pool.close()
                parallel_pool.join()

            else:  # for linux
                jobs = []
                queue = multiprocessing.Queue()
                worker_count = 1
                for lipid_sub in lipid_sub_key_lst:
                    if isinstance(lipid_sub, tuple) or isinstance(lipid_sub, list):

                        lipid_sub_lst = lipid_sub[0]
                        _chk_info_df = lipid_sub[1]
                        _chk_info_gp = lipid_sub[2]
                        if None in lipid_sub_lst:
                            lipid_sub_lst = [x for x in lipid_sub_lst if x is not None]
                        else:
                            pass
                        if isinstance(lipid_sub_lst[0], tuple) or isinstance(lipid_sub_lst[0], list):
                            lipid_sub_dct = {k: lipid_spec_dct[k] for k in lipid_sub_lst}
                        else:
                            lipid_sub_dct = {lipid_sub_lst: lipid_spec_dct[lipid_sub_lst]}
                            lipid_sub_lst = tuple([lipid_sub_lst])
                        print('[STATUS] >>> Core #%i ==> ...... processing ......' % worker_count)
                        job = multiprocessing.Process(target=get_lipid_info, args=(param_dct, usr_fa_df,
                                                                                   _chk_info_df,
                                                                                   _chk_info_gp, lipid_sub_lst,
                                                                                   usr_weight_df, key_frag_dct,
                                                                                   lipid_sub_dct, xic_dct,
                                                                                   worker_count,
                                                                                   save_fig, os_typ, queue))

                        worker_count += 1
                        jobs.append(job)
                        job.start()
                        lipid_info_results_lst.append(queue.get())

                for j in jobs:
                    j.join()

            # Merge multiprocessing results
            for lipid_info_result in lipid_info_results_lst:
                # TODO (georgia.angelidou@uni-leipzig.de):
                #  need to change the following part. The differnt sections can be combine

                if os_typ == 'windows':
                    try:
                        tmp_lipid_info = lipid_info_result.get()
                        # TODO (georgia.angelidou@uni-leipzig.de): can be avoided
                        tmp_lipid_info_df = tmp_lipid_info[0]
                        tmp_lipid_img_lst = tmp_lipid_info[1]

                        # TODO (georgia.angelidou@uni-leipzig.de): when above remove this need to be activated
                        # if isinstance(tmp_lipid_info[0], pd.DataFrame):
                        #     if not  tmp_lipid_info[0].empty:
                        #         output_df = output_df.append(tmp_lipid_info[0])
                        #         lipid_info_img_lst.extend(tmp_lipid_info[1])
                    except (KeyError, SystemError, ValueError, TypeError):
                        tmp_lipid_info_df = 'error'
                        tmp_lipid_img_lst = []
                        print('[ERROR] !!! This segment receive no Lipid identified.')
                else:  # for linux
                    try:
                        tmp_lipid_info_df = lipid_info_result[0]
                        tmp_lipid_img_lst = lipid_info_result[1]
                        # TODO (georgia.angelidou@uni-leipzig.de):
                        #  need to be activated when the above 2 lines are removed
                        # if isinstance(lipid_info_result[0], pd.DataFrame):
                        #     if not lipid_info_result[0].empty:
                        #         output_df = output_df.append(lipid_info_result[0])
                        #         lipid_info_img_lst.extend(lipid_info_result[1])

                    except (KeyError, SystemError, ValueError, TypeError):
                        tmp_lipid_info_df = 'error'
                        tmp_lipid_img_lst = []
                        print('[ERROR] !!! This segment receive no Lipid identified.')

                # TODO (georgia.angelidou@uni-leipzig.de):
                #  when new section are activate all the below if/else statement cant go away
                if isinstance(tmp_lipid_info_df, str):
                    pass
                else:
                    if isinstance(tmp_lipid_info_df, pd.DataFrame):
                        if not tmp_lipid_info_df.empty:
                            output_df = output_df.append(tmp_lipid_info_df)
                            lipid_info_img_lst.extend(tmp_lipid_img_lst)
            if part_tot == 1:
                print('[STATUS] >>> multiprocessing results merged')
            else:
                print('[STATUS] >>> multiprocessing results merged ==> Part %i / %i '
                      % (part_counter, part_tot))
            part_counter += 1

    else:
        print('[INFO] --> Using single core mode...')
        for lipid_sub_key_lst in lipid_part_key_lst:
            for lipid_sub in lipid_sub_key_lst:
                if isinstance(lipid_sub, tuple) or isinstance(lipid_sub, list):

                    lipid_sub_lst = lipid_sub[0]
                    _chk_info_df = lipid_sub[1]
                    _chk_info_gp = lipid_sub[2]

                    worker_count = 1
                    lipid_info_results_lst = get_lipid_info(param_dct, usr_fa_df, _chk_info_df,
                                                            _chk_info_gp, found_spec_key_lst, usr_weight_df,
                                                            key_frag_dct, lipid_spec_dct, xic_dct, worker_count)

                    # TODO (georgia.angelidou@uni-leipzig.de): unnecessary declaration of values can be avoided
                    tmp_lipid_info_df = lipid_info_results_lst[0]
                    tmp_lipid_img_lst = lipid_info_results_lst[1]
                    if isinstance(tmp_lipid_info_df, pd.DataFrame):
                        if not tmp_lipid_info_df.empty:
                            output_df = output_df.append(tmp_lipid_info_df)
                            lipid_info_img_lst = tmp_lipid_img_lst
                    # TODO (georgia.angelidou@uni-leipzig.de): code can change as below

    print('[OUTPUT] ==> Generate the output table')
    if isinstance(output_df, pd.DataFrame):
        print('[INFO] --> Total number of records', output_df.shape[0])
    if not output_df.empty:
        try:
            output_df = output_df.sort_values(by=['Lib_mz', 'Bulk_identification', 'MS2_scan_time', 'RANK_SCORE'],
                                              ascending=[True, True, True, False])
        except KeyError:
            pass
        output_df.reset_index(drop=True, inplace=True)
        output_df.index += 1

        output_df.drop_duplicates(keep='first', inplace=True)
        output_header_lst = output_df.columns.values.tolist()
        # TODO (georgia.angeldou@uni-leipzig.de): Add the info for the DG
        if usr_lipid_class in ['PA', 'PC', 'PE', 'PG', 'PI', 'PIP', 'PS']:
            output_list = ['FA1_[FA-H]-_i', 'FA2_[FA-H]-_i', '[LPL(FA1)-H]-_i', '[LPL(FA2)-H]-_i',
                           '[LPL(FA1)-H2O-H]-_i', '[LPL(FA2)-H2O-H]-_i']
            output_round_dct = {r'MS1_obs_mz': 4, r'Lib_mz': 4, 'ppm': 2, 'MS2_scan_time': 3,
                                'i_fa1': 2, 'i_fa2': 2, 'i_[M-H]-fa1': 2, 'i_[M-H]-fa2': 2,
                                'i_[M-H]-fa1-H2O': 2, 'i_[M-H]-fa2-H2O': 2
                                }

        elif usr_lipid_class in ['LPA', 'LPC', 'LPE', 'LPG', 'LPI', 'LPIP', 'LPS']:
            output_list = ['FA1_[FA-H]-_i']
            output_round_dct = {r'MS1_obs_mz': 4, r'Lib_mz': 4, 'ppm': 2, 'MS2_scan_time': 3,
                                'i_fa1': 2, 'i_fa2': 2, 'i_[M-H]-fa1': 2, 'i_[M-H]-fa2': 2,
                                'i_[M-H]-fa1-H2O': 2, 'i_[M-H]-fa2-H2O': 2
                                }

        elif usr_lipid_class in ['TG'] and usr_charge in ['[M+NH4]+', '[M+H]+']:
            output_list = ['FA1_[FA-H2O+H]+_i', 'FA2_[FA-H2O+H]+_i', 'FA3_[FA-H2O+H]+_i', '[MG(FA1)-H2O+H]+_i',
                           '[MG(FA2)-H2O+H]+_i', '[MG(FA3)-H2O+H]+_i', '[M-(FA1)+H]+_i', '[M-(FA2)+H]+_i',
                           '[M-(FA3)+H]+_i']
            output_round_dct = {r'MS1_obs_mz': 4, r'Lib_mz': 4, 'ppm': 2, 'MS2_scan_time': 3,
                                'i_fa1': 2, 'i_fa2': 2, 'i_fa3': 2, 'i_[M+H]-fa1': 2, 'i_[M+H]-fa2': 2,
                                'i_[M+H]-fa3': 2, 'i_[MG(fa1)+H]-H2O': 2, 'i_[MG(fa2)+H]-H2O': 2,
                                'i_[MG(fa3)+H]-H2O': 2}
        elif usr_lipid_class in ['TG'] and usr_charge in ['[M+Na]+']:
            # TODO (georgia.angelidou@uni-leipzig.de): check why this can cause some problems with [MGSN1-H2O+H]+
            output_list = ['FA1_[FA-H2O+H]+_i', 'FA2_[FA-H2O+H]+_i', 'FA3_[FA-H2O+H]+_i', '[MG(FA1)-H2O+H]+_i',
                           '[MG(FA2)-H2O+H]+_i', '[MG(FA3)-H2O+H]+_i', '[M-(FA1)+Na]+_i', '[M-(FA2)+Na]+_i',
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

        if len(target_ident_lst) > 0:
            for _t in target_ident_lst:
                output_round_dct[_t] = 2
        output_df = output_df.round(output_round_dct)

        output_df.rename(columns={'OBS_RESIDUES': '#Observed_FA'},
                         inplace=True)
        if usr_lipid_class in ['PA', 'PC', 'PE', 'PG', 'PI', 'PIP', 'PS']:
            output_short_lst = ['Proposed_structures', 'DISCRETE_ABBR', 'Formula_neutral', 'Formula_ion', 'Charge',
                                'Lib_mz', 'ppm', 'ISOTOPE_SCORE', 'RANK_SCORE',
                                'MS1_obs_mz', 'MS1_obs_i', r'MS2_PR_mz', 'MS2_scan_time',
                                'DDA#', 'Scan#', '#Observed_FA', '#Specific_peaks', '#Unspecific_peaks',
                                'FA1_[FA-H]-_i', 'FA2_[FA-H]-_i',
                                '[LPL(FA1)-H]-_i', '[LPL(FA2)-H]-_i',
                                '[LPL(FA1)-H2O-H]-_i', '[LPL(FA2)-H2O-H]-_i'
                                ]
        elif usr_lipid_class in ['LPA', 'LPC', 'LPE', 'LPG', 'LPI', 'LPIP', 'LPS']:
            output_short_lst = ['Proposed_structures', 'DISCRETE_ABBR', 'Formula_neutral', 'Formula_ion', 'Charge',
                                'Lib_mz', 'ppm', 'ISOTOPE_SCORE', 'RANK_SCORE',
                                'MS1_obs_mz', 'MS1_obs_i', r'MS2_PR_mz', 'MS2_scan_time',
                                'DDA#', 'Scan#', '#Observed_FA', '#Specific_peaks', '#Unspecific_peaks',
                                'FA1_[FA-H]-_i',
                                ]
        elif usr_lipid_class in ['TG'] and usr_charge in ['[M+H]+', '[M+NH4]+']:
            output_short_lst = ['Proposed_structures', 'DISCRETE_ABBR', 'Formula_neutral', 'Formula_ion', 'Charge',
                                'Lib_mz', 'ppm', 'ISOTOPE_SCORE', 'RANK_SCORE',
                                'MS1_obs_mz', 'MS1_obs_i', r'MS2_PR_mz', 'MS2_scan_time',
                                'DDA#', 'Scan#', '#Observed_FA',
                                '[M-(FA1)+H]+_i', '[M-(FA2)+H]+_i', '[M-(FA3)+H]+_i',
                                '[MG(FA1)-H2O+H]+_i', '[MG(FA2)-H2O+H]+_i', '[MG(FA3)-H2O+H]+_i',
                                'FA1_[FA-H2O+H]+_i', 'FA2_[FA-H2O+H]+_i', 'FA3_[FA-H2O+H]+_i',
                                ]
        elif usr_lipid_class in ['TG'] and usr_charge in ['[M+Na]+']:
            output_short_lst = ['Proposed_structures', 'DISCRETE_ABBR', 'Formula_neutral', 'Formula_ion', 'Charge',
                                'Lib_mz', 'ppm', 'ISOTOPE_SCORE', 'RANK_SCORE',
                                'MS1_obs_mz', 'MS1_obs_i', r'MS2_PR_mz', 'MS2_scan_time',
                                'DDA#', 'Scan#', '#Observed_FA',
                                '[M-(FA1)+Na]+_i', '[M-(FA2)+Na]+_i', '[M-(FA3)+Na]+_i',
                                '[MG(FA1)-H2O+H]+_i', '[MG(FA2)-H2O+H]+_i', '[MG(FA3)-H2O+H]+_i',
                                'FA1_[FA-H2O+H]+_i', 'FA2_[FA-H2O+H]+_i', 'FA3_[FA-H2O+H]+_i'
                                ]
        elif usr_lipid_class in ['DG'] and usr_charge in ['[M+H]+', '[M+NH4]+', '[M+Na]+']:
            output_short_lst = ['Proposed_structures', 'DISCRETE_ABBR', 'Formula_neutral', 'Formula_ion', 'Charge',
                                'Lib_mz', 'ppm', 'ISOTOPE_SCORE', 'RANK_SCORE',
                                'MS1_obs_mz', 'MS1_obs_i', r'MS2_PR_mz', 'MS2_scan_time',
                                'DDA#', 'Scan#', '#Observed_FA',
                                '[MG(FA1)-H2O+H]+_i', '[MG(FA2)-H2O+H]+_i',
                                'FA1_[FA-H2O+H]+_i', 'FA2_[FA-H2O+H]+_i',
                                ]
        else:
            output_short_lst = ['Proposed_structures', 'DISCRETE_ABBR', 'Formula_neutral', 'Formula_ion', 'Charge',
                                'Lib_mz', 'ppm', 'ISOTOPE_SCORE', 'RANK_SCORE',
                                'MS1_obs_mz', 'MS1_obs_i', r'MS2_PR_mz', 'MS2_scan_time',
                                'DDA#', 'Scan#', '#Observed_FA',
                                ]

        # Make copy of selected list of columns from output_df to avoid pandas warnings
        final_output_df = output_df[output_short_lst]
        final_output_df = final_output_df.sort_values(by=['MS1_obs_mz', 'MS2_scan_time', 'RANK_SCORE'],
                                                      ascending=[True, True, False])
        final_output_df = final_output_df.reset_index(drop=True)
        final_output_df.index += 1
        save_output(output_sum_xlsx, final_output_df)

    else:
        error_lst.append('[Warning] NO Lipid identified in this file.\n!! Please check your settings !!')
        tot_run_time = time.clock() - start_time
        print('[WARNING] !!! This file got no Lipid identified.')
        print('[STATUS] >>> Identification finished in %s sec <<<' % tot_run_time)
        return tot_run_time, error_lst, output_df

    # Reserved for the development to re-use the session
    if save_session is True:
        hunt_save_path = os.path.join(output_folder, 'HunterData_%s.hunt' % hunter_start_time_str)
        results_pickle_dct = {'param_dct': param_dct, 'output_df': output_df,
                              'final_output_df': final_output_df, 'lipid_info_img_lst': lipid_info_img_lst}
        save_hunt(results_pickle_dct, hunt_save_path)
        print('Hunter session saved as:', hunt_save_path)
    else:
        pass

    # Start multiprocessing to save img for HTML report
    if save_fig is True:
        gen_html_report(param_dct, output_df, lipid_info_img_lst)
    else:
        print('[WARNING] !!! User skip image generation !!!!!!')

    tot_run_time = time.clock() - start_time

    print('[STATUS] >>> >>> >>> FINISHED in %s sec <<< <<< <<<' % tot_run_time)

    return tot_run_time, error_lst, output_df
