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

import pickle

import math
import multiprocessing
from multiprocessing import Pool
import os
import time

try:
    from LibLipidHunter.LogPageCreator import LogPageCreator
    from LibLipidHunter.PanelPlotter import gen_plot
except ImportError:  # for python 2.7.14
    from LogPageCreator import LogPageCreator
    from PanelPlotter import gen_plot


def save_hunt(results_pickle_dct, hunt_save_path):

    for key in ['lipid_info_img_lst', 'param_dct', 'output_df', 'final_output_df']:
        if key in list(results_pickle_dct.keys()):
            pass
        else:
            print('Key %s is missing' % key)
            raise KeyError

    with open(hunt_save_path, 'wb') as results_pickle:
        pickle.dump(results_pickle_dct, results_pickle)


def recover_hunt(hunter_data):

    hunt_pickle_dct = pickle.load(hunter_data)

    print(list(hunt_pickle_dct.keys()))
    param_dct = hunt_pickle_dct['param_dct']
    output_df = hunt_pickle_dct['output_df']
    lipid_info_img_lst = hunt_pickle_dct['lipid_info_img_lst']
    final_output_df = hunt_pickle_dct['final_output_df']

    usr_vendor = param_dct['vendor']
    output_folder = param_dct['img_output_folder_str']
    output_sum_xlsx = param_dct['xlsx_output_path_str']
    usr_ms1_ppm = param_dct['ms_ppm']
    usr_ms2_ppm = param_dct['ms2_ppm']
    usr_ms1_precision = usr_ms1_ppm * 1e-6
    hunter_start_time_str = param_dct['hunter_start_time']
    usr_core_num = param_dct['core_number']
    usr_dpi = param_dct['img_dpi']
    usr_img_type = param_dct['img_type']

    output_sum_xlsx_directory = os.path.dirname(output_sum_xlsx)
    if not os.path.exists(output_sum_xlsx_directory):
        os.makedirs(output_sum_xlsx_directory)
    try:
        final_output_df.to_excel(output_sum_xlsx, index=False)
        print('[OUTPUT] ==> Prepare to save output as: ', output_sum_xlsx)
    except IOError:
        final_output_df.to_excel('%s-%i%s' % (output_sum_xlsx[:-5], int(time.time()), '.xlsx'), index=False)
        print(output_sum_xlsx)
    print('[OUTPUT] ==> File saved ...')

    current_path = os.getcwd()
    if os.path.isdir(output_folder):
        os.chdir(output_folder)
        if os.path.isdir('LipidHunter_Results_Figures_%s' % hunter_start_time_str):
            print('[INFO] --> Output folder existed...')
        else:
            os.mkdir('LipidHunter_Results_Figures_%s' % hunter_start_time_str)
            print('[INFO] --> Output folder created...')
    else:
        os.mkdir(output_folder)
        os.chdir(output_folder)
        os.mkdir('LipidHunter_Results_Figures_%s' % hunter_start_time_str)
        print('[INFO] --> Output folder created...')
    os.chdir(current_path)

    # generate html files
    log_pager = LogPageCreator(output_folder, hunter_start_time_str, param_dct)
    log_pager.add_all_info(output_df)
    log_pager.close_page()
    # del log_pager
    print('[STATUS] >>> start to generate images: image count %i' % len(lipid_info_img_lst))

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

                if len(img_sub_lst) > 0:
                    print('[STATUS] >>> Core #%i ==> Generating output images ... image count: %i'
                          % (worker_count, len(img_sub_lst)))
                    if 'debug_mode' in list(param_dct.keys()):
                        if param_dct['debug_mode'] == 'ON':
                            # TODO (georgia.angelidou@uni-leipzig.de): Check if the following is really necesary or can be done with alternative way
                            for img_param_dct in img_sub_lst:
                                print(img_param_dct['save_img_as'])
                            del img_param_dct
                    parallel_pool.apply_async(gen_plot, args=(img_sub_lst, worker_count, usr_img_type,
                                                              usr_dpi, usr_vendor, usr_ms1_precision))
                    worker_count += 1

        parallel_pool.close()
        parallel_pool.join()

    else:
        worker_count = 1
        print('[INFO] --> Using single core mode...')
        if isinstance(lipid_info_img_lst, tuple) or isinstance(lipid_info_img_lst, list):
            if None in lipid_info_img_lst:
                lipid_info_img_lst = [x for x in lipid_info_img_lst if x is not None]
            else:
                pass
            if len(lipid_info_img_lst) > 0:
                gen_plot(lipid_info_img_lst, worker_count, usr_img_type, usr_dpi,
                         usr_vendor, usr_ms1_precision)


if __name__ == '__main__':

    from test.test_HuntManager import test_recover_hunt

    test_recover_hunt()
