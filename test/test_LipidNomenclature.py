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

from LibLipidHunter.LipidNomenclature import NameParserFA

cwd = os.getcwd()
if cwd.endswith('test'):
    print('change to folder above..')
    os.chdir('..')

print('test_get_fa_info @ ', os.getcwd())


def test_get_fa_info():
    abbr_decoder = NameParserFA()
    abbr_lst = ['FA16:0', 'FA18:0', 'FA18:1', 'O-16:0', 'P-18:0']
    for abbr in abbr_lst:
        x = abbr_decoder.get_fa_info(abbr)
        print(x)
