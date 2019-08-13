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

import logging
import os
import sys
import unittest

hunterPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, hunterPath + "/../")

from LibLipidHunter.LipidComposer import LipidComposer

log_level = logging.DEBUG
logging.basicConfig(
    format="%(asctime)s-%(levelname)s - %(message)s",
    datefmt="%b-%d@%H:%M:%S",
    level=log_level,
)
logger = logging.getLogger("log")


class TestCaseLipidComposer(unittest.TestCase):
    def setUp(self):
        logger.debug("SETUP TESTS... TestCase_cmd_lipidhunter")
        cwd = os.getcwd()
        if cwd.endswith("test") or cwd.endswith("test/") or cwd.endswith("test\\"):
            logger.info("change to folder above..")
            os.chdir("..")
        logger.info(os.getcwd())

        print("TestCaseLipidComposer @ ", os.getcwd())
        self.fa_lst_file = r"ConfigurationFiles/1-FA_Whitelist.xlsx"

        # Note:
        # exact position means to consider the position from the FA white list that the user give but,
        # in the case that the user define 2 different FA for both positions then:
        # When it is false it will give only one option
        # and when it is TRUE to give both combinations that these 2 FA an make (in case of phospholipids)

        self.bad_param_dct = {
            "fa_whitelist": r"file_not_exist_FA_list.xlsx",
            "lipid_class": "lx",
            "charge_mode": "[M+x]-",
            "exact_position": "FALSE",
        }
        self.lpc_param_dct = {
            "fa_whitelist": self.fa_lst_file,
            "lipid_class": "LPC",
            "charge_mode": "[M+HCOO]-",
            "exact_position": "FALSE",
        }
        self.lpe_param_dct = {
            "fa_whitelist": self.fa_lst_file,
            "lipid_class": "LPE",
            "charge_mode": "[M-H]-",
            "exact_position": "FALSE",
        }
        self.pc_param_dct = {
            "fa_whitelist": self.fa_lst_file,
            "lipid_class": "PC",
            "charge_mode": "[M+HCOO]-",
            "exact_position": "FALSE",
        }
        self.pe_param_dct = {
            "fa_whitelist": self.fa_lst_file,
            "lipid_class": "PE",
            "charge_mode": "[M-H]-",
            "exact_position": "FALSE",
        }
        self.dg_param_dct = {
            "fa_whitelist": self.fa_lst_file,
            "lipid_class": "DG",
            "charge_mode": "[M+NH4]+",
            "exact_position": "FALSE",
        }
        self.tg_param_dct = {
            "fa_whitelist": self.fa_lst_file,
            "lipid_class": "TG",
            "charge_mode": "[M+NH4]+",
            "exact_position": "FALSE",
        }
        # Todo: add SM parameters here
        if not os.path.isdir(r"test/results/LipidComposer"):
            os.makedirs(r"test/results/LipidComposer")

    def get_lipidmaster(self, lipid_dct):
        is_successful = False
        logger.info(f"Test LipidComposer for {lipid_dct['lipid_class']}")
        composer = LipidComposer()
        lipid_master_df = composer.compose_lipid(param_dct=lipid_dct, ms2_ppm=30)
        master_csv = (
            r"test/results/LipidComposer/LipidMaster_Whitelist_%s.csv"
            % lipid_dct["lipid_class"]
        )
        fa_csv = (
            r"test/results/LipidComposer/LipidMaster_FAlist_%s.csv"
            % lipid_dct["lipid_class"]
        )
        try:
            os.remove(master_csv)
            os.remove(fa_csv)
            logger.info("Previous left over data deleted ...")
        except Exception as e:
            logger.info(e)
            logger.info("No previous left over data ...")
        calc_fa_df = composer.calc_fa_query(
            lipid_dct["lipid_class"], self.fa_lst_file, ms2_ppm=50
        )
        if not calc_fa_df.empty and not lipid_master_df.empty:
            logger.debug(calc_fa_df.head())
            lipid_master_df.to_csv(master_csv)
            calc_fa_df.to_csv(fa_csv)
            if os.path.isfile(master_csv) and os.path.isfile(fa_csv):
                is_successful = True

        return is_successful

    def test_bad_file(self):
        is_sucessful = True
        try:
            is_sucessful = self.get_lipidmaster(self.bad_param_dct)
        except Exception as e:
            logger.error(e)
            is_sucessful = False
        assert is_sucessful is False

    def test_lpc(self):
        is_sucessful = self.get_lipidmaster(self.lpc_param_dct)
        assert is_sucessful is True

    def test_lpe(self):
        is_sucessful = self.get_lipidmaster(self.lpe_param_dct)
        assert is_sucessful is True

    def test_pc(self):
        is_sucessful = self.get_lipidmaster(self.pc_param_dct)
        assert is_sucessful is True

    def test_pe(self):
        is_sucessful = self.get_lipidmaster(self.pe_param_dct)
        assert is_sucessful is True

    # Todo: add SM test function here

    def test_dg(self):
        is_sucessful = self.get_lipidmaster(self.tg_param_dct)
        assert is_sucessful is True

    def test_tg(self):
        is_sucessful = self.get_lipidmaster(self.tg_param_dct)
        assert is_sucessful is True

    def tearDown(self):
        logger.debug("Test LipidComposer TEST PASSED!")


if __name__ == "__main__":
    unittest.main()
    logger.info("TESTS FINISHED!")
