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

import logging
import os
import sys
import unittest

hunterPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, hunterPath + '/../')

import subprocess
import cmd_lipidhunter

log_level = logging.DEBUG
logging.basicConfig(format='%(asctime)s-%(levelname)s - %(message)s', datefmt='%b-%d@%H:%M:%S', level=log_level)
logger = logging.getLogger('log')

import pytest


class TestCase_cmd_lipidhunter(unittest.TestCase):

    def setUp(self):
        logger.debug('SETUP TESTS... TestCase_cmd_lipidhunter')
        cwd = os.getcwd()
        if cwd.endswith('test') or cwd.endswith('test/') or cwd.endswith('test\\'):
            logger.info('change to folder above..')
            os.chdir('..')
        logger.info(os.getcwd())
        pl_cfg_path = r'test/test_batch_cfg/test_PC_cfg.txt'
        tg_cfg_path = r'test/test_batch_cfg/test_TG_cfg.txt'
        bad_cfg_path = r'badtest/test_batch_cfg/test_bad_cfg.txt'

        self.pass_params_pl = ['-i', pl_cfg_path]
        self.pass_params_tg = ['-i', tg_cfg_path]
        self.fail_input_params = ['-i', bad_cfg_path]

    def test_cmd_lipidhunter_help(self):
        logger.debug('Test help...')
        assert cmd_lipidhunter.main(['-h']) is False

    def test_cmd_lipidhunter_bad_params(self):
        logger.debug('Test bad params...')
        assert cmd_lipidhunter.main(['-test']) is False

    def test_cmd_lipidhunter_bad_infile(self):
        logger.debug('Test bad input...')
        assert cmd_lipidhunter.main(self.fail_input_params) is False

    def test_cmd_lipidhunter_pl(self):
        logger.debug('Test sample data... PL')
        assert cmd_lipidhunter.main(self.pass_params_pl) is True

    # @pytest.mark.skip(reason="Skip data recover for Travis CI")
    def test_cmd_lipidhunter_tg(self):
        logger.debug('Test sample data ... TG')
        assert cmd_lipidhunter.main(self.pass_params_tg) is True

    def tearDown(self):
        logger.debug('TestCase_cmd_lipidhunter TEST END!')


if __name__ == '__main__':
    # python cmd_lipidhunter.py -i test/test_batch_cfg/test_PC_cfg.txt
    unittest.main()
    logger.info('TESTS FINISHED!')
