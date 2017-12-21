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
import sys

import pandas as pd
from multiprocessing import Pool

from ParallelFunc import ppm_calc_para, ppm_window_para, pr_window_calc_para


def find_pr_info(scan_info_df, spectra_pl, lpp_info_groups, sub_group_list, ms1_th, ms1_ppm, ms1_max):
    core_results_df = pd.DataFrame()
    for group_key in sub_group_list:
        subgroup_df = lpp_info_groups.get_group(group_key)

        same_mz_se = subgroup_df.iloc[0, :].squeeze()
        _pr_code = '%f<= MS2_PR_mz <= %f' % (same_mz_se['PR_MZ_LOW'], same_mz_se['PR_MZ_HIGH'])

        _tmp_scan_info_df = scan_info_df.query(_pr_code).copy()

        if _tmp_scan_info_df.shape[0] > 0:

            for idx, row in _tmp_scan_info_df.iterrows():
                ms2_pr = row['MS2_PR_mz']
                ms2_dda_idx = row['dda_event_idx']
                ms2_dda_rank = row['DDA_rank']
                ms2_spec_idx = row['spec_index']
                ms2_scan_time = row['scan_time']
                ms2_scan_number = row['scan_number']

                _ms1_query_code = 'dda_event_idx == %i and DDA_rank == 0' % ms2_dda_idx
                tmp_ms1_info_df = scan_info_df.query(_ms1_query_code).copy()

                if tmp_ms1_info_df.shape[0] > 0:

                    ms1_spec_idx = tmp_ms1_info_df['spec_index'].tolist()[0]

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

                        if pr_ms1_df.shape[0] > 0:

                            pr_ms1_df.loc[:, 'ppm'] = ppm_calc_para(pr_ms1_df['mz'].values, same_mz_se['Lib_mz'])
                            pr_ms1_df.loc[:, 'abs_ppm'] = abs(pr_ms1_df.loc[:, 'ppm'])
                            pr_info_df = pr_ms1_df.query('abs_ppm <= %i' % ms1_ppm).copy()
                            if pr_info_df.shape[0] > 0:
                                pr_info_df.sort_values(by='i', ascending=False, inplace=True)
                                pr_info_df.reset_index(drop=True, inplace=True)
                                _ms1_mz = pr_info_df['mz'].tolist()[0]
                                _ppm = pr_info_df['ppm'].tolist()[0]
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
                                print('core_results_df.shape', core_results_df.shape)
                            else:
                                print('pr_info_df.shape[0] == 0')
                        else:
                            print('pr_ms1_df.shape[0] == 0')
                    else:
                        print('ms1_spec_idx not in list')
                else:
                    print('tmp_ms1_info_df.shape[0] = 0')
        else:
            print('_tmp_scan_info_df.shape[0] = 0')

    print(core_results_df.shape)

    return core_results_df


class PrecursorHunter(object):
    def __init__(self, lpp_info_df, param_dct):
        self.lpp_info_df = lpp_info_df
        self.param_dct = param_dct

    def get_matched_pr(self, scan_info_df, spectra_pl, ms1_max=0, core_num=4, max_ram=8):

        print('Start match precursors!!!!')

        pr_window = self.param_dct['pr_window']

        ms1_ppm = self.param_dct['ms_ppm']
        ms1_th = self.param_dct['ms_th']

        pl_class = self.param_dct['lipid_type']
        charge_mode = self.param_dct['charge_mode']

        ms1_obs_pr_df = pd.DataFrame()

        if pl_class == 'PC' and charge_mode == '[M+HCOO]-':

            pc_fa_df = self.lpp_info_df[self.lpp_info_df['[M+HCOO]-_MZ'] > 0]
            pc_h_df = self.lpp_info_df[self.lpp_info_df['[M-H]-_MZ'] > 0]

            print('PC [M+HCOO]- LPP: ', pc_fa_df.shape[0])
            print('PC [M-H]- LPP: ', pc_h_df.shape[0])
            pc_fa_df.loc[:, 'PR_MZ_LOW'] = pr_window_calc_para(pc_fa_df['[M+HCOO]-_MZ'].values.tolist(), -1 * pr_window)
            pc_fa_df.loc[:, 'PR_MZ_HIGH'] = pr_window_calc_para(pc_fa_df['[M+HCOO]-_MZ'].values.tolist(), pr_window)
            pc_fa_df.loc[:, 'MS1_MZ_LOW'] = ppm_window_para(pc_fa_df['[M+HCOO]-_MZ'].values.tolist(), -1 * ms1_ppm)
            pc_fa_df.loc[:, 'MS1_MZ_HIGH'] = ppm_window_para(pc_fa_df['[M+HCOO]-_MZ'].values.tolist(), ms1_ppm)
            pc_fa_df.loc[:, 'Formula'] = pc_fa_df['[M+HCOO]-_FORMULA'].str.strip('-')
            pc_fa_df.loc[:, 'Ion'] = '[M+HCOO]-'
            pc_fa_df.loc[:, 'Lib_mz'] = pc_fa_df.loc[:, '[M+HCOO]-_MZ']
            if pc_h_df.shape[0] > 0:
                pc_h_df.loc[:, 'PR_MZ_LOW'] = pr_window_calc_para(pc_h_df['[M-H]-_MZ'].values.tolist(), -1 * pr_window)
                pc_h_df.loc[:, 'PR_MZ_HIGH'] = pr_window_calc_para(pc_h_df['[M-H]-_MZ'].values.tolist(), pr_window)
                pc_h_df.loc[:, 'MS1_MZ_LOW'] = ppm_window_para(pc_h_df['[M-H]-_MZ'].values.tolist(), -1 * ms1_ppm)
                pc_h_df.loc[:, 'MS1_MZ_HIGH'] = ppm_window_para(pc_h_df['[M-H]-_MZ'].values.tolist(), ms1_ppm)
                pc_h_df.loc[:, 'Formula'] = pc_h_df['[M-H]-_FORMULA'].str.strip('-')
                pc_h_df.loc[:, 'Ion'] = '[M-H]-'
                pc_h_df.loc[:, 'Lib_mz'] = pc_h_df.loc[:, '[M-H]-_MZ']
                self.lpp_info_df = pc_fa_df.copy()
                self.lpp_info_df = self.lpp_info_df.append(pc_h_df)
            else:
                self.lpp_info_df = pc_fa_df.copy()

            self.lpp_info_df = self.lpp_info_df.sort_values(by='PR_MZ_LOW')

        else:
            self.lpp_info_df.loc[:, 'PR_MZ_LOW'] = pr_window_calc_para(self.lpp_info_df['[M-H]-_MZ'].values.tolist(),
                                                                       -1 * pr_window)
            self.lpp_info_df.loc[:, 'PR_MZ_HIGH'] = pr_window_calc_para(self.lpp_info_df['[M-H]-_MZ'].values.tolist(),
                                                                        pr_window)
            self.lpp_info_df.loc[:, 'MS1_MZ_LOW'] = ppm_window_para(self.lpp_info_df['[M-H]-_MZ'].values.tolist(),
                                                                    -1 * ms1_ppm)
            self.lpp_info_df.loc[:, 'MS1_MZ_HIGH'] = ppm_window_para(self.lpp_info_df['[M-H]-_MZ'].values.tolist(),
                                                                     ms1_ppm)
            self.lpp_info_df.loc[:, 'Formula'] = self.lpp_info_df.loc[:, '[M-H]-_FORMULA'].str.strip('-')
            self.lpp_info_df.loc[:, 'Ion'] = '[M-H]-'
            self.lpp_info_df.loc[:, 'Lib_mz'] = self.lpp_info_df.loc[:, '[M-H]-_MZ']

        # Prepare for multiprocessing
        lpp_info_groups = self.lpp_info_df.groupby(['Lib_mz', 'Formula'])
        all_group_key_lst = lpp_info_groups.groups.keys()
        sub_len = int(math.ceil(len(all_group_key_lst) / core_num))
        core_key_list = map(None, *(iter(all_group_key_lst),) * sub_len)

        spectra_pl_idx_lst = sorted(spectra_pl.items.tolist())

        if len(spectra_pl_idx_lst) >= (max_ram * 50):
            sub_pl_group_lst = map(None, *(iter(spectra_pl_idx_lst),) * (max_ram * 40))
        else:
            sub_pl_group_lst = [spectra_pl_idx_lst]

        # print('sub_pl_group_lst')
        # print(sub_pl_group_lst)

        part_tot = len(sub_pl_group_lst)
        part_counter = 1
        opt_sub_pl_group_lst = []
        for sub_idx_lst in sub_pl_group_lst:

            if isinstance(sub_idx_lst, tuple) or isinstance(sub_idx_lst, list):
                sub_idx_lst = filter(lambda x: x is not None, sub_idx_lst)
                opt_sub_pl_group_lst.append(sub_idx_lst)
                sub_pl = spectra_pl.loc[sub_idx_lst, :, :]
                print(sub_pl.items)

                # Start multiprocessing
                if part_tot == 1:
                    print('>>> Start multiprocessing for precursor matching ==> Number of Cores: %i' % core_num)
                else:
                    print('>>> Start multiprocessing for precursor matching ==> Part %i / %i --> Number of Cores: %i' %
                          (part_counter, part_tot, core_num))

                if self.param_dct['core_number'] > 1:
                    part_counter += 1
                    parallel_pool = Pool(core_num)
                    pr_info_results_lst = []
                    core_worker_count = 1
                    for core_list in core_key_list:
                        if isinstance(core_list, tuple) or isinstance(core_list, list):
                            if None in core_list:
                                core_list = filter(lambda x: x is not None, core_list)
                            else:
                                pass
                            print('>>> >>> ...... Core #%i ==> processing ......' % core_worker_count)
                            pr_info_result = parallel_pool.apply_async(find_pr_info, args=(scan_info_df,
                                                                                           sub_pl,
                                                                                           lpp_info_groups,
                                                                                           core_list, ms1_th,
                                                                                           ms1_ppm, ms1_max))
                            core_worker_count += 1
                            pr_info_results_lst.append(pr_info_result)

                    parallel_pool.close()
                    parallel_pool.join()

                    for pr_info_result in pr_info_results_lst:
                        try:
                            sub_df = pr_info_result.get()
                            if sub_df.shape[0] > 0:
                                ms1_obs_pr_df = ms1_obs_pr_df.append(sub_df)
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
                            print('>>> >>> processing ......Part: %i subset: %i ' % (part_counter, core_worker_count))
                            sub_df = find_pr_info(scan_info_df, sub_pl, lpp_info_groups, core_list, ms1_th, ms1_ppm, ms1_max)
                            core_worker_count += 1
                            if sub_df.shape[0] > 0:
                                ms1_obs_pr_df = ms1_obs_pr_df.append(sub_df)
                    part_counter += 1

        # End multiprocessing

        print('ms1_obs_pr_df.shape', ms1_obs_pr_df.shape)
        if ms1_obs_pr_df.shape[0] > 0:
            ms1_obs_pr_df = ms1_obs_pr_df.sort_values(by=['Lib_mz', 'abs_ppm'], ascending=[True, True])
            ms1_obs_pr_df = ms1_obs_pr_df.reset_index(drop=True)
            return ms1_obs_pr_df, opt_sub_pl_group_lst

        else:
            return '!! NO suitable precursor --> Check settings!!', 'Empty spec_lst'

