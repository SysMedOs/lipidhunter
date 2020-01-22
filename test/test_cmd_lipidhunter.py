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
import subprocess
import unittest

import pytest

hunterPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, hunterPath + "/../")

import cmd_lipidhunter

log_level = logging.DEBUG
logging.basicConfig(
    format="%(asctime)s-%(levelname)s - %(message)s",
    datefmt="%b-%d@%H:%M:%S",
    level=log_level,
)
logger = logging.getLogger("log")


class TestCase_cmd_lipidhunter(unittest.TestCase):
    def setUp(self):
        logger.debug("SETUP TESTS... TestCase_cmd_lipidhunter")
        cwd = os.getcwd()
        if cwd.endswith("test") or cwd.endswith("test/") or cwd.endswith("test\\"):
            logger.info("change to folder above..")
            os.chdir("..")
        logger.info(os.getcwd())
        if not os.path.isdir(r"test/results"):
            os.makedirs(r"test/results")
            logger.info(f"folder created: {os.getcwd()}/test/results")
        lpc_cfg_path = r"test/test_batch_cfg/test_LPC_cfg.txt"
        lpe_cfg_path = r"test/test_batch_cfg/test_LPE_cfg.txt"
        pc_cfg_path = r"test/test_batch_cfg/test_PC_cfg.txt"
        pe_cfg_path = r"test/test_batch_cfg/test_PE_cfg.txt"
        dg_cfg_path = r"test/test_batch_cfg/test_DG_cfg.txt"
        tg_NH4_cfg_path = r"test/test_batch_cfg/test_TG_NH4_cfg.txt"
        tg_Na_cfg_path = r"test/test_batch_cfg/test_TG_Na_cfg.txt"
        tg_H_cfg_path = r"test/test_batch_cfg/test_TG_H_cfg.txt"
        bad_cfg_path = r"badtest/test_batch_cfg/test_bad_cfg.txt"

        self.pass_params_lpc = ["-i", lpc_cfg_path]
        self.pass_params_lpe = ["-i", lpe_cfg_path]
        self.pass_params_pc = ["-i", pc_cfg_path]
        self.pass_params_pe = ["-i", pe_cfg_path]
        self.pass_params_dg = ["-i", dg_cfg_path]
        self.pass_params_tg_NH4 = ["-i", tg_NH4_cfg_path]
        self.pass_params_tg_Na = ["-i", tg_Na_cfg_path]
        self.pass_params_tg_H = ["-i", tg_H_cfg_path]
        self.fail_input_params = ["-i", bad_cfg_path]

    def test_help(self):
        logger.debug("Test help...")
        assert cmd_lipidhunter.main(["-h"]) is False

    def test_bad_params(self):
        logger.debug("Test bad params...")
        assert cmd_lipidhunter.main(["-test"]) is False

    def test_bad_infile(self):
        logger.debug("Test bad input...")
        assert cmd_lipidhunter.main(self.fail_input_params) is False

    @pytest.mark.skip(reason="Currently no Lyso PL files")
    def test_lpc(self):
        logger.debug("Test sample data... PL")
        assert cmd_lipidhunter.main(self.pass_params_lpc) is True

    @pytest.mark.skip(reason="Currently no Lyso PL files")
    def test_lpe(self):
        logger.debug("Test sample data... PL")
        assert cmd_lipidhunter.main(self.pass_params_lpe) is True

    def test_pc(self):
        logger.debug("Test sample data... PL")
        assert cmd_lipidhunter.main(self.pass_params_pc) is True

    def test_pe(self):
        logger.debug("Test sample data... PL")
        assert cmd_lipidhunter.main(self.pass_params_pe) is True

    @pytest.mark.skip(reason="Currently no DG files")
    def test_dg(self):
        logger.debug("Test sample data ... DG")
        assert cmd_lipidhunter.main(self.pass_params_dg) is True

    def test_tg_NH4(self):
        logger.debug("Test sample data ... TG [M+NH4]+")
        assert cmd_lipidhunter.main(self.pass_params_tg_NH4) is True

    @pytest.mark.skip(reason="Currently no TG Na files")
    def test_tg_Na(self):
        logger.debug("Test sample data ... TG [M+Na]+")
        assert cmd_lipidhunter.main(self.pass_params_tg_Na) is True

    @pytest.mark.skip(reason="Currently no TG H files")
    def test_tg_H(self):
        logger.debug("Test sample data ... TG [M+H]+")
        assert cmd_lipidhunter.main(self.pass_params_tg_H) is True

    def tearDown(self):
        logger.debug("TestCase_cmd_lipidhunter TEST PASSED!")


if __name__ == "__main__":
    # python cmd_lipidhunter.py -i test/test_batch_cfg/test_PC_cfg.txt
    unittest.main()
    logger.info("TESTS FINISHED!")
