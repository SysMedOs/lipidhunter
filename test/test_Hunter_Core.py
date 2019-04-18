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

import os
import sys
import time
import unittest

hunterPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, hunterPath + '/../')

from LibLipidHunter.Hunter_Core import huntlipids

cwd = os.getcwd()
if cwd.endswith('test'):
    print('change to folder above..')
    os.chdir('..')


def test_huntlipids():
    # set the core number and max ram in GB to be used for the test
    core_count = 4
    max_ram = 5  # int only
    save_images = True  # True --> generate images, False --> NO images (not recommended)
    save_lipidmaster_table = True  # True --> export LipidMasterTable to output folder, False --> NO export
    save_hunter_session = True  # True --> export session as a file, False --> NO export

    # full_test_lst = [['PC', 'waters'],['PE', 'waters'], ['TG', 'waters','[M+H]+'], ['TG', 'waters', '[M+NH4]+'],
    # ['TG', 'thermo', '[M+NH4]+']]

    test_lst = [
        ['PC', 'waters', '[M+HCOO]-', 'PC_waters'],
        # ['PE', 'waters', '[M-H]-', 'PE_waters'],
        # ['TG', 'waters', '[M+H]+', 'TG_waters'],
        # ['TG', 'waters', '[M+NH4]+', 'TG_waters_NH4'],
        # ['TG', 'waters', '[M+Na]+', 'TG_waters_Na'],
        ['TG', 'thermo', '[M+NH4]+', 'TG_thermo_NH4'],
        # ['DG', 'thermo', '[M+NH4]+', 'DG_thermo_NH4'],
        # ['DG', 'agilent', '[M+NH4]+', 'TG_agilent_NH4'],
        # ['LPC', 'thermo', '[M+HCOO]-', 'LPC_thermo'],
        # ['LPE', 'thermo', '[M-H]-', 'LPE_thermo'],
        # ['LPA', 'thermo', '[M-H]-', 'LPA_thermo'],
        # ['LPG', 'thermo', '[M-H]-', 'LPG_thermo'],
        # ['LPS', 'thermo', '[M-H]-', 'LPS_thermo'],
        # ['LPI', 'thermo', '[M-H]-', 'LPI_thermo'],
    ]

    # set the default files
    pl_mzml_waters = r'test/mzML/PL_Neg_Waters_qTOF.mzML'  # Synapt-g2si file
    lpl_mzml_thermo = r'test/mzML/dev_test/LPL_Neg_Thermo_Orbi.mzML'  # Qexactive file
    tg_mzml_waters = r'test/mzML/dev_test/TG_Pos_Waters_qTOF.mzML'  # Synapt-g2si file
    tg_mzml_thermo = r'test/mzML/TG_Pos_Thermo_Orbi.mzML'  # Qexactive file
    tg_mzml_SCIEXS = r'test/mzML/Test_sciex.mzML'  # position holder
    tg_mzml_agilent = r'test/mzML/Test_agilent.mzML'  # position holder

    pl_base_dct = {'fawhitelist_path_str': r'ConfigurationFiles/1-FA_Whitelist.xlsx',
                   'lipid_specific_cfg': r'ConfigurationFiles/3-Specific_ions.xlsx',
                   'score_cfg': r'ConfigurationFiles/2-Score_weight_PL.xlsx'}

    lpl_base_dct = {'fawhitelist_path_str': r'ConfigurationFiles/1-FA_Whitelist.xlsx',
                    'lipid_specific_cfg': r'ConfigurationFiles/3-Specific_ions.xlsx',
                    'score_cfg': r'ConfigurationFiles/2-Score_weight_LPL.xlsx'}

    tg_base_dct = {'fawhitelist_path_str': r'ConfigurationFiles/1-FA_Whitelist.xlsx',
                   'lipid_specific_cfg': r'ConfigurationFiles/3-Specific_ions.xlsx',
                   'score_cfg': r'ConfigurationFiles/2-Score_weight_TG.xlsx'}

    usr_test_dct = {}
    usr_test_dct_keys = []
    mz_range = [600, 1000]  # default
    rt_range = [6, 10]  # default
    for usr_test in test_lst:
        _test_dct = {'rank_score_filter': 60.5, 'score_filter': 60.5, 'isotope_score_filter': 75.0, 'ms_max': 0,
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
                mz_range = [600, 1000]  # 600, 1000
                rt_range = [24, 27]  # max [24, 27]
            else:
                mzml = False
                pass

            _test_dct.update(pl_base_dct)

        elif usr_test[0] in ['LPC', 'LPE', 'LPA', 'LPG', 'LPI', 'LPS', 'LPIP']:
            lipid_class = usr_test[0]
            if lipid_class == 'LPC':
                charge = '[M+HCOO]-'
            else:
                charge = '[M-H]-'
            vendor = usr_test[1]
            if vendor == 'waters':
                mzml = pl_mzml_waters
                mz_range = [600, 1000]  # 600, 1000
                rt_range = [24, 27]  # max [24, 27]
            elif vendor == 'thermo':
                mzml = lpl_mzml_thermo
                mz_range = [300, 900]  # 600, 1000 short 800 - 840
                rt_range = [3, 10]  # [20, 28] short [24, 26]
            else:
                mzml = False
                pass

            _test_dct.update(lpl_base_dct)

        elif usr_test[0] in ['TG']:
            lipid_class = usr_test[0]
            vendor = usr_test[1]
            charge = usr_test[2]
            if vendor == 'waters':
                mzml = tg_mzml_waters
                mz_range = [800, 1000]  # 600, 1000
                rt_range = [9, 15]  # max [9, 15]
            elif vendor == 'thermo':
                mzml = tg_mzml_thermo
                mz_range = [600, 1000]  # 600, 1000 short 800 - 840
                rt_range = [20, 28]  # [20, 28] short [24, 26]
            elif vendor == 'sciex':
                mzml = tg_mzml_SCIEXS
                mz_range = [900, 1000]  # 600, 1000
                rt_range = [8, 13]
            elif vendor == 'agilent':
                mzml = tg_mzml_agilent
                mz_range = [600, 1000]  # 600, 1000
                rt_range = [9, 13]
            else:
                mzml = False

            if lipid_class == 'TG' and charge == '[M+Na]+':
                tg_base_dct['score_cfg'] = r'ConfigurationFiles/2-Score_weight_TG_Na.xlsx'

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
                mz_range = [300, 850]
                rt_range = [5, 15]
            elif vendor == 'agilent':
                mzml = tg_mzml_agilent
                mz_range = [400, 1000]
                rt_range = [5, 13]
            else:
                mzml = False

            _test_dct.update(tg_base_dct)
        else:

            lipid_class = False
            charge = False
            vendor = False
            mzml = False

        if mzml is not False and charge is not False:

            _cfg_dct = {'lipid_class': lipid_class, 'charge_mode': charge, 'vendor': vendor, 'mzml_path_str': mzml,
                        'rt_start': rt_range[0], 'rt_end': rt_range[1], 'mz_start': mz_range[0], 'mz_end': mz_range[1]}

            if vendor == 'waters':
                _cfg_dct['ms_ppm'] = 20
                _cfg_dct['ms2_ppm'] = 50  # 50
                _cfg_dct['ms_th'] = 750
                _cfg_dct['ms2_th'] = 10
                _cfg_dct['dda_top'] = 6

                _test_dct.update(_cfg_dct)
                usr_test_dct[usr_test[3]] = _test_dct
                usr_test_dct_keys.append(usr_test[3])

            elif vendor == 'thermo':
                _cfg_dct['ms_ppm'] = 10
                _cfg_dct['ms2_ppm'] = 50
                _cfg_dct['ms_th'] = 500
                _cfg_dct['ms2_th'] = 50
                _cfg_dct['dda_top'] = 30  # 10

                _test_dct.update(_cfg_dct)
                usr_test_dct[usr_test[3]] = _test_dct
                usr_test_dct_keys.append(usr_test[3])
            elif vendor == 'sciex':
                _cfg_dct['ms_ppm'] = 10
                _cfg_dct['ms2_ppm'] = 60
                _cfg_dct['ms_th'] = 5000
                _cfg_dct['ms2_th'] = 100  # Can get 2000/1000/750/500 depends how strict should be the identification
                _cfg_dct['dda_top'] = 15
                ms_ppm_SCIEXS = 10

                _test_dct.update(_cfg_dct)
                usr_test_dct[usr_test[3]] = _test_dct
                usr_test_dct_keys.append(usr_test[3])
            elif vendor == 'agilent':
                _cfg_dct['ms_ppm'] = 50
                _cfg_dct['ms2_ppm'] = 100
                _cfg_dct['ms_th'] = 1000
                _cfg_dct['ms2_th'] = 10  # Can get 2000/1000/750/500 depends how strict should be the identification
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
    hunter_folder = os.getcwd()
    hunter_file_path = os.path.join(hunter_folder, 'LipidHunter.py')

    if os.path.isfile(hunter_file_path):
        print('\nLipidHunter folder', hunter_folder, '\n')
        print(usr_test_dct, '\n')

        for test_key in usr_test_dct_keys:
            if test_key in list(usr_test_dct.keys()):
                test_dct = usr_test_dct[test_key]
                t_str = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
                lipid_class = test_dct['lipid_class']

                cfg_dct = {'img_output_folder_str': r'Test/results/%s_%s' % (test_key, t_str),
                           'xlsx_output_path_str': r'Test/results/%s_%s.xlsx' % (test_key, t_str),
                           'hunter_folder': hunter_folder, 'img_type': u'png', 'img_dpi': 300,
                           'hunter_start_time': t_str, 'experiment_mode': 'LC-MS', 'rank_score': True,
                           'fast_isotope': False, 'core_number': core_count, 'max_ram': max_ram, 'tag_all_sn': True}

                test_dct.update(cfg_dct)
                if save_lipidmaster_table is True:
                    test_dct['debug_mode'] = 'ON'
                    test_dct['save_lipid_master_table'] = 'CSV'
                print(test_dct)

                print('>>>>>>>>>>>>>>>> START TEST: %s' % test_key)

                t, log_lst, output_df2 = huntlipids(test_dct, log_lst,
                                                    save_fig=save_images, save_session=save_hunter_session)
                if t is not False:
                    print('>>>>>>>>>>>>>>>> TEST PASSED: %s in %.3f Sec <<<<<<<<<<<<<<<<\n' % (test_key, t))
                    t_sum_lst.append((test_key, 'PASSED', '%.3f Sec' % t, 'Identified:'))
                else:
                    print('>>>>>>>>!!!!!!!! TEST FAILED: %s !!!!!!!<<<<<<<<\n' % test_key)
                    t_sum_lst.append((test_key, 'FAILED', '', 'Identified: 0'))

    else:
        print('!!! Invalid LipidHunter folder', hunter_folder)
        print('Please run HunterCore from LibLipidHunter folder')
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
