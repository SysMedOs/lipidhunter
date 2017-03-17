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
#
# For more info please contact:
#     SysMedOs_team: oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
#     Developer Georgia Angelidou georgia.angelidou@uni-leipzig.de
#


import os
import sys
import multiprocessing
from PySide import QtGui
from LibLipidHunter.BatchHunter_Core import BatchHunterCore


if __name__ == '__main__':
    multiprocessing.freeze_support()
    usr_cwd = os.getcwd()
    gui = QtGui.QApplication(sys.argv)
    BatchLipidHunter = BatchHunterCore(cwd=usr_cwd)
    BatchLipidHunter.show()
