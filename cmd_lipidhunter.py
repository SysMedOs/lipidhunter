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

from __future__ import print_function

import ConfigParser as configparser
import getopt
import os.path
import sys
import time

from LibLipidHunter.Hunter_Core import huntlipids


def main(argv):
    """
    To run LipidHunter from command line, please generate one configuration file by GUI mode and use it as a template.
    You can load each time one configuration file only.
    :param argv: -i <input LipidHunter configuration file in .txt format>
    """

    _cfg_file = ''
    cfg_params_dct = {}

    i_type_key_lst = ['ms_th', 'ms2_th', 'hg_th', 'ms_ppm', 'ms2_ppm', 'hg_ppm', 'dda_top', 'sn_ratio',
                      'core_number', 'max_ram', 'img_dpi', 'ms_max']
    f_type_key_lst = ['rt_start', 'rt_end', 'mz_start', 'mz_end', 'pr_window', 'ms2_infopeak_threshold',
                      'ms2_hginfopeak_threshold', 'score_filter', 'isotope_score_filter', 'rank_score_filter']
    b_type_key_lst = ['rank_score', 'fast_isotope', 'tag_all_sn']

    try:
        opts, args = getopt.getopt(argv, 'hi:o:', ['infile='])
    except getopt.GetoptError:
        print('Error: cmd_lipidhunter.py -i <input LipidHunter configuration file in .txt format>')
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            print('cmd_lipidhunter.py -i <input LipidHunter configuration file in .txt format>')
            sys.exit(1)
        elif opt in ('-i', '--infile'):
            _cfg_file = arg

    if isinstance(_cfg_file, str) and len(_cfg_file) > 0:
        print('Input LipidHunter configuration file : ', _cfg_file)
        if os.path.isfile(_cfg_file):
            with open(_cfg_file) as _cfg_obj:
                config = configparser.ConfigParser()
                config.readfp(_cfg_obj)
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
                    sys.exit(1)
        else:
            print('Load configuration file !!! FAILED !!! File do not exist !!!')
            sys.exit(1)
    else:
        print('Error: cmd_lipidhunter.py -i <input LipidHunter configuration file in .txt format>')
        sys.exit(1)

    if len(cfg_params_dct.keys()) > 0:
        start_time_str = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
        cfg_params_dct['hunter_start_time'] = start_time_str
        print('Sart to process... ', start_time_str)
        t, log_lst = huntlipids(cfg_params_dct, error_lst=[])
        if len(log_lst) > 0:
            for err in log_lst:
                print('Error:', err)
            print('Error: Please check your parameters!!!')
        else:
            print(t)
            print('test passed!')


if __name__ == "__main__":
    main(sys.argv[1:])
