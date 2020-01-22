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

from LibLipidHunter.LipidNomenclature import NameParserFA

log_level = logging.DEBUG
logging.basicConfig(
    format="%(asctime)s-%(levelname)s - %(message)s",
    datefmt="%b-%d@%H:%M:%S",
    level=log_level,
)
logger = logging.getLogger("log")


class TestCase_LipidNomenclature(unittest.TestCase):
    def setUp(self):
        logger.debug("SETUP TESTS... TestCase_LipidNomenclature")
        cwd = os.getcwd()
        if cwd.endswith("test") or cwd.endswith("test/") or cwd.endswith("test\\"):
            logger.info("change to folder above..")
            os.chdir("..")
        logger.info(os.getcwd())
        self.pass_params = ["FA16:0", "FA18:0", "FA18:1", "O-16:0", "P-18:0"]
        self.bad_params = ["BAD", "guys", "here", "x4!"]

    def test_bad_params(self):
        logger.debug("Test bad params...")
        abbr_decoder = NameParserFA()
        result_lst = []
        for abbr in self.bad_params:
            try:
                x = abbr_decoder.get_fa_info(abbr)
                print(x)
                if x:
                    result_lst.append(x)
            except Exception as e:
                logger.error(f"Can not parse {abbr}")
                logger.error(f"Error {e}")
        assert len(result_lst) < len(self.bad_params)

    def test_good_params(self):
        logger.debug("Test good params...")
        abbr_decoder = NameParserFA()
        result_lst = []
        for abbr in self.pass_params:
            x = abbr_decoder.get_fa_info(abbr)
            print(x)
            if x:
                result_lst.append(x)
        assert len(result_lst) == len(self.pass_params)

    def tearDown(self):
        logger.debug("TestCase_cmd_lipidhunter TEST PASSED!")


if __name__ == "__main__":
    # python cmd_lipidhunter.py -i test/test_batch_cfg/test_PC_cfg.txt
    unittest.main()
    logger.info("TESTS FINISHED!")
