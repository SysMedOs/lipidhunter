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

import math
import multiprocessing
from multiprocessing import Pool
from sys import platform

import pandas as pd

try:
    from LibLipidHunter.ParallelFunc import ppm_calc_para, ppm_window_para, pr_window_calc_para

except ImportError:  # for python 2.7
    from ParallelFunc import ppm_calc_para, ppm_window_para, pr_window_calc_para


def find_pr_info(scan_info_df, spectra_pl, lpp_info_groups, sub_group_list, ms1_th, ms1_ppm, ms1_max, core=1,
                 os_type='windows', queue=None):
    core_count = 'Core #{core}'.format(core=core)
    print(core_count, '... Matching precursors ...')
    core_results_df = pd.DataFrame()
    for group_key in sub_group_list:
        subgroup_df = lpp_info_groups.get_group(group_key).copy()
        subgroup_df.is_copy = False

        same_mz_se = subgroup_df.iloc[0, :].squeeze()
        _pr_code = '%f<= MS2_PR_mz <= %f' % (same_mz_se['PR_MZ_LOW'], same_mz_se['PR_MZ_HIGH'])

        _tmp_scan_info_df = scan_info_df.query(_pr_code).copy()

        if not _tmp_scan_info_df.empty:

            for idx, row in _tmp_scan_info_df.iterrows():
                ms2_pr = row['MS2_PR_mz']
                ms2_dda_idx = row['dda_event_idx']
                ms2_dda_rank = row['DDA_rank']
                ms2_spec_idx = row['spec_index']
                ms2_scan_time = row['scan_time']
                ms2_scan_number = row['scan_number']
                _ms1_query_code = 'dda_event_idx == %i and DDA_rank == 0' % ms2_dda_idx
                tmp_ms1_info_df = scan_info_df.query(_ms1_query_code).copy()

                if not tmp_ms1_info_df.empty:

                    ms1_spec_idx = tmp_ms1_info_df['spec_index'].values.tolist()[0]

                    if ms1_spec_idx in spectra_pl.items:
                        ms1_df = spectra_pl[ms1_spec_idx]
                        if ms1_max > ms1_th:
                            pr_ms1_df = ms1_df.query('%f <= i <= %f and %f <= mz <= %f' % (ms1_th, ms1_max,
                                                                                           same_mz_se['MS1_MZ_LOW'],
                                                                                           same_mz_se['MS1_MZ_HIGH']))
                        else:
                            pr_ms1_df = ms1_df.query('%f <= i and %f <= mz <= %f' % (ms1_th,
                                                                                     same_mz_se['MS1_MZ_LOW'],
                                                                                     same_mz_se['MS1_MZ_HIGH']))

                        pr_ms1_df.is_copy = False

                        if not pr_ms1_df.empty:

                            pr_ms1_df.loc[:, 'ppm'] = ppm_calc_para(pr_ms1_df['mz'].values, same_mz_se['Lib_mz'])
                            pr_ms1_df.loc[:, 'abs_ppm'] = abs(pr_ms1_df.loc[:, 'ppm'])
                            pr_info_df = pr_ms1_df.query('abs_ppm <= %i' % ms1_ppm).copy()
                            if not pr_info_df.empty:
                                pr_info_df.sort_values(by='i', ascending=False, inplace=True)
                                pr_info_df.reset_index(drop=True, inplace=True)
                                _ms1_mz = pr_info_df['mz'].values.tolist()[0]
                                _ppm = pr_info_df['ppm'].values.tolist()[0]
                                len_df = subgroup_df.shape[0]
                                subgroup_df.loc[:, 'MS1_obs_mz'] = [_ms1_mz] * len_df
                                subgroup_df.loc[:, 'dda_event_idx'] = [ms2_dda_idx] * len_df
                                subgroup_df.loc[:, 'spec_index'] = [ms2_spec_idx] * len_df
                                subgroup_df.loc[:, 'scan_time'] = [ms2_scan_time] * len_df
                                subgroup_df.loc[:, 'DDA_rank'] = [ms2_dda_rank] * len_df
                                subgroup_df.loc[:, 'scan_number'] = [ms2_scan_number] * len_df
                                subgroup_df.loc[:, 'MS2_PR_mz'] = [ms2_pr] * len_df
                                subgroup_df.loc[:, 'MS1_XIC_mz'] = [round(_ms1_mz, 4)] * len_df
                                subgroup_df.loc[:, 'ppm'] = [_ppm] * len_df
                                subgroup_df.loc[:, 'abs_ppm'] = [abs(_ppm)] * len_df

                                core_results_df = core_results_df.append(subgroup_df)
                                # print('core_results_df.shape', core_results_df.shape)
                            else:
                                pass
                                # print('pr_info_df.shape[0] == 0')
                        else:
                            pass
                            # print('pr_ms1_df.shape[0] == 0')
                    else:
                        pass
                        # print('ms1_spec_idx not in list')
                else:
                    pass
                    # print('tmp_ms1_info_df.shape[0] = 0')
        else:
            pass
            # print('_tmp_scan_info_df.shape[0] = 0')

    print(core_count, 'core_results_count', core_results_df.shape[0])

    if os_type == 'linux_multi':
        queue.put(core_results_df)
    else:
        return core_results_df


class PrecursorHunter(object):
    def __init__(self, lpp_info_df, param_dct, os_type='windows'):
        self.lpp_info_df = lpp_info_df.copy()
        self.lpp_info_df.is_copy = False
        self.param_dct = param_dct
        self.os_typ = os_type

    def get_matched_pr(self, scan_info_df, spectra_pl, ms1_max=0, core_num=4, max_ram=8):

        print('>>>>>>>> Start match precursors!!!!')

        if core_num > 4:
            core_num = 4
            print('>>> Temporarily reduce to 4 cores to enhance the performance in this step...')

        pr_window = self.param_dct['pr_window']

        ms1_ppm = self.param_dct['ms_ppm']
        ms1_th = self.param_dct['ms_th']

        pl_class = self.param_dct['lipid_type']
        usr_charge = self.param_dct['charge_mode']

        ms1_obs_pr_df = pd.DataFrame()
        # print('ms1_obs_pr_df')
        if not self.lpp_info_df.empty:
            if pl_class == 'PC':
                if usr_charge in ['[M+HCOO]-', '[M+CH3COO]-']:
                    lpp_mz_lst = self.lpp_info_df['%s_MZ' % usr_charge].values.tolist()
                    self.lpp_info_df.loc[:, 'PR_MZ_LOW'] = pr_window_calc_para(lpp_mz_lst, -1 * pr_window)
                    self.lpp_info_df.loc[:, 'PR_MZ_HIGH'] = pr_window_calc_para(lpp_mz_lst, pr_window)
                    self.lpp_info_df.loc[:, 'MS1_MZ_LOW'] = ppm_window_para(lpp_mz_lst, -1 * ms1_ppm)
                    self.lpp_info_df.loc[:, 'MS1_MZ_HIGH'] = ppm_window_para(lpp_mz_lst, ms1_ppm)
                    self.lpp_info_df.loc[:, 'Formula'] = self.lpp_info_df['%s_FORMULA' % usr_charge].str.strip('-')
                    self.lpp_info_df.loc[:, 'Ion'] = usr_charge
                    self.lpp_info_df.loc[:, 'Lib_mz'] = self.lpp_info_df.loc[:, '%s_MZ' % usr_charge]
                else:
                    self.lpp_info_df = pd.DataFrame()

            else:
                lpp_mz_lst = self.lpp_info_df['%s_MZ' % usr_charge].values.tolist()
                self.lpp_info_df.loc[:, 'PR_MZ_LOW'] = pr_window_calc_para(lpp_mz_lst, -1 * pr_window)
                self.lpp_info_df.loc[:, 'PR_MZ_HIGH'] = pr_window_calc_para(lpp_mz_lst, pr_window)
                self.lpp_info_df.loc[:, 'MS1_MZ_LOW'] = ppm_window_para(lpp_mz_lst, -1 * ms1_ppm)
                self.lpp_info_df.loc[:, 'MS1_MZ_HIGH'] = ppm_window_para(lpp_mz_lst, ms1_ppm)
                self.lpp_info_df.loc[:, 'Formula'] = self.lpp_info_df.loc[:, '%s_FORMULA' % usr_charge].str.strip('-')
                self.lpp_info_df.loc[:, 'Ion'] = usr_charge
                self.lpp_info_df.loc[:, 'Lib_mz'] = self.lpp_info_df.loc[:, '%s_MZ' % usr_charge]

            self.lpp_info_df.sort_values(by='PR_MZ_LOW', inplace=True)

        else:
            # no matched info --> exit
            return False, False

        # Prepare for multiprocessing
        lpp_info_groups = self.lpp_info_df.groupby(['Lib_mz', 'Formula'])
        all_group_key_lst = list(lpp_info_groups.groups.keys())
        sub_len = int(math.ceil(len(all_group_key_lst) / core_num))
        core_key_list = [all_group_key_lst[k: k + sub_len] for k in range(0, len(all_group_key_lst), sub_len)]

        spectra_pl_idx_lst = sorted(spectra_pl.items.values.tolist())

        if (max_ram * 64) <= len(spectra_pl_idx_lst) < (max_ram * 128):
            print('>>>>>>>> Spectra is too large for the RAM settings, split to few segments ...')
            sub_group_len = int(math.ceil(len(spectra_pl_idx_lst) * 0.5))
        elif (max_ram * 128) <= len(spectra_pl_idx_lst) < (max_ram * 256):
            print('>>>>>>>> Spectra is too large for the RAM settings, split to few segments ...')
            sub_group_len = int(math.ceil(len(spectra_pl_idx_lst) * 0.25))
        elif len(spectra_pl_idx_lst) >= (max_ram * 256):
            print('>>>>>>>> Spectra is too large for the RAM settings, split to few segments ...')
            sub_group_len = int(math.ceil(len(spectra_pl_idx_lst) * 0.1))

        else:
            sub_group_len = len(spectra_pl_idx_lst)

        print('Total Scans:', len(spectra_pl_idx_lst), 'Sub part scans:', sub_group_len)
        if sub_group_len > 500:
            print('>>> Set sub part scans to 500 to avoid Memory error ...')
            sub_group_len = 500
        sub_pl_group_lst = [spectra_pl_idx_lst[s: (s + sub_group_len)] for s in range(0, len(spectra_pl_idx_lst),
                                                                                      sub_group_len)]

        part_tot = len(sub_pl_group_lst)
        part_counter = 1
        opt_sub_pl_group_lst = []
        pr_info_results_lst = []
        for sub_idx_lst in sub_pl_group_lst:

            if isinstance(sub_idx_lst, tuple) or isinstance(sub_idx_lst, list):
                sub_idx_lst = [x for x in sub_idx_lst if x is not None]
                opt_sub_pl_group_lst.append(sub_idx_lst)
                sub_pl = spectra_pl.loc[sub_idx_lst, :, :]
                # print(sub_pl.items)

                # Start multiprocessing
                if part_tot == 1:
                    print('>>> Start multiprocessing for precursor matching ==> Number of Cores: %i' % core_num)
                else:
                    print('>>> Start multiprocessing for precursor matching ==> Part %i / %i --> Number of Cores: %i' %
                          (part_counter, part_tot, core_num))

                if self.param_dct['core_number'] > 1:
                    if self.os_typ == 'windows':
                        parallel_pool = Pool(core_num)

                        core_worker_count = 1
                        for core_list in core_key_list:
                            if isinstance(core_list, tuple) or isinstance(core_list, list):
                                if None in core_list:
                                    core_list = [x for x in core_list if x is not None]
                                else:
                                    pass
                                print('>>> >>> ...... Core #%i ==> processing ......' % core_worker_count)
                                pr_info_result = parallel_pool.apply_async(find_pr_info, args=(scan_info_df,
                                                                                               sub_pl,
                                                                                               lpp_info_groups,
                                                                                               core_list, ms1_th,
                                                                                               ms1_ppm, ms1_max,
                                                                                               core_worker_count,
                                                                                               self.os_typ))
                                core_worker_count += 1
                                pr_info_results_lst.append(pr_info_result)

                        parallel_pool.close()
                        parallel_pool.join()
                    else:
                        jobs = []
                        queue = multiprocessing.Queue()
                        core_worker_count = 1
                        for core_list in core_key_list:
                            if isinstance(core_list, tuple) or isinstance(core_list, list):
                                if None in core_list:
                                    core_list = [x for x in core_list if x is not None]
                                else:
                                    pass
                                print('>>> >>> ...... Core #%i ==> processing ......' % core_worker_count)
                                job = multiprocessing.Process(target=find_pr_info, args=(scan_info_df, sub_pl,
                                                                                         lpp_info_groups, core_list,
                                                                                         ms1_th, ms1_ppm, ms1_max,
                                                                                         core_worker_count,
                                                                                         self.os_typ, queue,))
                                core_worker_count += 1
                                jobs.append(job)
                                job.start()
                                pr_info_results_lst.append(queue.get())
                        for j in jobs:
                            j.join()

                else:
                    print('Using single core mode...')
                    core_worker_count = 1
                    for core_list in core_key_list:
                        if isinstance(core_list, tuple) or isinstance(core_list, list):
                            if None in core_list:
                                core_list = [x for x in core_list if x is not None]
                            else:
                                pass
                            print('>>> >>> processing ......Part: %i subset: %i ' % (part_counter, core_worker_count))
                            sub_df = find_pr_info(scan_info_df, sub_pl, lpp_info_groups, core_list, ms1_th,
                                                  ms1_ppm, ms1_max, core_worker_count)
                            if not sub_df.empty:
                                pr_info_results_lst.append(sub_df)

                part_counter += 1

        #  Merge multiprocessing results
        result_counter = 0
        result_part_counter = 1
        for pr_info_result in pr_info_results_lst:
            if self.param_dct['core_number'] > 1:
                if self.os_typ == 'windows':
                    try:
                        sub_df = pr_info_result.get()
                        if not sub_df.empty:
                            ms1_obs_pr_df = ms1_obs_pr_df.append(sub_df)
                    except (KeyError, SystemError, ValueError, TypeError):
                        pass
                else:
                    try:
                        if not pr_info_result.empty:
                            ms1_obs_pr_df = ms1_obs_pr_df.append(pr_info_result)
                    except (KeyError, SystemError, ValueError, TypeError):
                        pass
            else:
                if not pr_info_result.empty:
                    ms1_obs_pr_df = ms1_obs_pr_df.append(pr_info_result)

            result_counter += 1
            if result_counter > core_num * result_part_counter:
                result_part_counter += 1

            if part_tot == 1:
                print('>>> Multiprocessing results merged ...')
            else:
                print('>>> Multiprocessing results merged ... Part %i / %i ...' % (result_part_counter, part_tot))

        # End multiprocessing

        print('ms1_obs_pr_df.shape', ms1_obs_pr_df.shape)
        if not ms1_obs_pr_df.empty:
            ms1_obs_pr_df = ms1_obs_pr_df.sort_values(by=['Lib_mz', 'abs_ppm'], ascending=[True, True])
            ms1_obs_pr_df = ms1_obs_pr_df.reset_index(drop=True)
            return ms1_obs_pr_df, opt_sub_pl_group_lst

        else:
            return False, False
