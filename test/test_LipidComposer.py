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
import unittest

hunterPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, hunterPath + '/../')
cwd = os.getcwd()
if cwd.endswith('test'):
    print('change to folder above..')
    os.chdir('..')

from LibLipidHunter.LipidComposer import LipidComposer

print('test_lipicomposer @ ', os.getcwd())


def test_lipicomposer():
    fa_lst_file = r'ConfigurationFiles/1-FA_Whitelist.xlsx'

    # Note:
    # exact position means to consider the poition from the FA white list that the user give but,
    # in the case that the user define 2 different FA for both positions then:
    # When it is false it will give only one option
    # and when it is TRUE to give both compinations that these 2 FA an make (in case of phospholipids)

    # usr_param_dct = {'fa_whitelist': fa_lst_file, 'lipid_class': 'TG', 'charge_mode': '[M+NH4]+',
    #                  'exact_position': 'FALSE'}

    # usr_param_dct = {'fa_whitelist': fa_lst_file, 'lipid_class': 'LPC', 'charge_mode': '[M+HCOO]-',
    #                  'exact_position': 'FALSE'}
    usr_param_dct = {'fa_whitelist': fa_lst_file, 'lipid_class': 'LPE', 'charge_mode': '[M-H]-',
                     'exact_position': 'FALSE'}

    composer = LipidComposer()
    usr_lipid_master_df = composer.compose_lipid(param_dct=usr_param_dct, ms2_ppm=30)
    print('[INFO] --> Lipid Master Table generated...')

    master_csv = r'Temp/LipidMaster_Whitelist_%s.csv' % usr_param_dct['lipid_class']
    fa_csv = r'Temp/LipidMaster_FAlist_%s.csv' % usr_param_dct['lipid_class']

    calc_fa_df = composer.calc_fa_query(usr_param_dct['lipid_class'], fa_lst_file, ms2_ppm=50)

    if calc_fa_df is False:
        print('[ERROR] !!! Failed to generate FA info table ...\n')

    print(calc_fa_df.head())
    usr_lipid_master_df.to_csv(master_csv)
    calc_fa_df.to_csv(fa_csv)
    print('[INFO] --> Finished...')
