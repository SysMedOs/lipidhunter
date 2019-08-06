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

import configparser
import getopt
# required to perform multiprocessing
import multiprocessing
import os.path
import sys
from sys import platform
import time
from typing import List

from LibLipidHunter.Hunter_Core import huntlipids


def main(argv: List[str]) -> bool:

    """
    To run LipidHunter from command line, please generate one configuration file by GUI mode and use it as a template.
    You can load each time one configuration file only.
    Args:
        argv (str): -i <input LipidHunter configuration file in .txt format>

    Returns:
        is_successful(bool): Return True if the run finished with no error message
    """

    is_successful = False

    _cfg_file = ''
    cfg_params_dct = {}

    i_type_key_lst = ['ms_th', 'ms2_th', 'hg_th', 'ms_ppm', 'ms2_ppm', 'hg_ppm', 'dda_top', 'sn_ratio',
                      'core_number', 'max_ram', 'img_dpi', 'ms_max']
    f_type_key_lst = ['rt_start', 'rt_end', 'mz_start', 'mz_end', 'pr_window', 'ms2_infopeak_threshold',
                      'ms2_hginfopeak_threshold', 'score_filter', 'isotope_score_filter', 'rank_score_filter']
    b_type_key_lst = ['rank_score', 'fast_isotope', 'tag_all_sn']

    save_img = True

    try:
        opts, args = getopt.getopt(argv, 'hi:o:n', ['infile='])
    except getopt.GetoptError:
        print('Error: cmd_lipidhunter.py -i <input LipidHunter configuration file in .txt format>')
        return is_successful

    for opt, arg in opts:
        if opt == '-h':
            print('python cmd_lipidhunter.py -i <input LipidHunter configuration file in .txt format>')
            print('Use -n to skip output image generation (not recommended).')
            return is_successful
        elif opt in ('-i', '--infile'):
            _cfg_file = arg
        elif opt == '-n':
            save_img = False

    if isinstance(_cfg_file, str) and len(_cfg_file) > 0:
        print('Input LipidHunter configuration file : ', _cfg_file)
        if os.path.isfile(_cfg_file):
            with open(_cfg_file) as _cfg_obj:
                config = configparser.ConfigParser()
                try:
                    config.read_file(_cfg_obj)
                except AttributeError:  # for python 2.7.14
                    config.readfp(_cfg_obj)
                print('got file', _cfg_file)
                if config.has_section('parameters'):
                    usr_cfg = 'parameters'
                    options = config.options(usr_cfg)
                    for param in options:
                        _val = config.get(usr_cfg, param)
                        if param in i_type_key_lst:
                            try:
                                cfg_params_dct[param] = int(_val)
                            except ValueError:
                                cfg_params_dct[param] = int(float(_val))
                        elif param in f_type_key_lst:
                            cfg_params_dct[param] = float(_val)
                        elif param in b_type_key_lst:
                            if _val.lower() == 'true':
                                cfg_params_dct[param] = True
                            if _val.lower() == 'false':
                                cfg_params_dct[param] = False
                        else:
                            cfg_params_dct[param] = _val
                    print('Load configuration file... Passed ...')
                else:
                    print('Error: Load configuration file FAILED !!! Configuration file content error !!!')
                    return is_successful
        else:
            print('Load configuration file !!! FAILED !!! File do not exist !!!')
            return is_successful
    else:
        print('Error: cmd_lipidhunter.py -i <input LipidHunter configuration file in .txt format>')
        return is_successful

    if len(list(cfg_params_dct.keys())) > 0:
        start_time_str = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
        cfg_params_dct['hunter_start_time'] = start_time_str

        output_folder_path = cfg_params_dct['img_output_folder_str']

        if platform == "linux" or platform == "linux2":
            l_cwd = os.getcwd()
            if output_folder_path.startswith('/'):
                os.chdir('/')
            if os.path.isdir(output_folder_path):
                print('Folder existed...\n', output_folder_path)
                folder_abs_path = os.path.abspath(output_folder_path)
                print('abs path of folder\n', folder_abs_path)
            else:
                if os.path.isdir('/' + output_folder_path):
                    print('Folder existed...\n', output_folder_path)
                    folder_abs_path = os.path.abspath(output_folder_path)
                    print('abs path of folder\n', folder_abs_path)
                else:
                    print('No folder...\n', output_folder_path)
                    os.makedirs(output_folder_path)
                    print('Folder created... %s' % output_folder_path)
            os.chdir(l_cwd)
        else:
            if os.path.isdir(output_folder_path):
                print('Output folder path... %s' % output_folder_path)
            else:
                try:
                    os.mkdir(output_folder_path)
                    print('Output folder created... %s' % output_folder_path)
                except IOError:
                    print('!! Failed to create output folder !!')
        print('Start to process... ', start_time_str)
        t, log_lst, export_df = huntlipids(cfg_params_dct, error_lst=[], save_fig=save_img)
        if len(log_lst) > 0:
            for err in log_lst:
                print('Error:', err)
            print('Error: Please check your parameters!!!')
        else:
            print(t)
            print('Run finished!')
            is_successful = True

    return is_successful


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main(sys.argv[1:])
