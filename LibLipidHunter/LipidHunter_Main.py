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

from __future__ import division
from __future__ import print_function

import glob
import multiprocessing
import multiprocessing.pool
import os
from sys import platform
import re
import time

import pandas as pd
from PySide import QtCore, QtGui
from six.moves import configparser

try:
    from LibLipidHunter.LipidHunter_UI import Ui_MainWindow
    from LibLipidHunter.Hunter_Core import huntlipids
    from LibLipidHunter.LipidComposer import LipidComposer
except ImportError:  # for python 2.7.14
    from LipidHunter_UI import Ui_MainWindow
    from Hunter_Core import huntlipids
    from LipidComposer import LipidComposer


class LipidHunterMain(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None, cwd=None):
        scale = 1
        config = configparser.ConfigParser()
        config.read('config.ini')
        if config.has_section('settings'):
            user_cfg = 'settings'
        else:
            if config.has_section('default'):
                user_cfg = 'default'
            else:
                user_cfg = ''
        if len(user_cfg) > 2:
            options = config.options(user_cfg)
            if 'gui_scale' in options:
                try:
                    scale = float(config.get(user_cfg, 'gui_scale'))
                except (ValueError, TypeError):
                    pass
        if scale != 1:
            if scale > 2.75:
                scale = 2.75
            elif scale < 0.75:
                scale = 0.75
            print('[INFO] Using GUI scale x{0}'.format(scale))
        QtGui.QMainWindow.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self, scale=scale)

        # set version
        version_date = r'06, March, 2019'
        version_html = (r'<html><head/><body><p><span style=" font-weight:600;">'
                        r'LipidHunter 2 (RC) # Released Date: {version_date}'
                        r'</span></p></body></html>').format(version_date=version_date)
        self.ui.version_lb.setText(QtGui.QApplication.translate("MainWindow", version_html, None,
                                                                QtGui.QApplication.UnicodeUTF8))

        # current folder:
        if cwd is not None:
            print('User LipidHunter folder', cwd)
            self.lipidhunter_cwd = cwd
        else:
            auto_cwd = os.getcwd()
            print('User LipidHunter folder', auto_cwd)
            self.lipidhunter_cwd = auto_cwd

        self.load_cfg()
        self.a_max_ms()

        # disable un necessary UI elements
        q_pgb_style = '''
                      QProgressBar:horizontal {
                      border: 1px #CCCCCC; border-radius: 0px; background: #CCCCCC; padding: 1px;}
                      QProgressBar::chunk:horizontal {
                      background: qlineargradient(x1: 0, y1: 0.5, x2: 0, y2: 0.5, stop: 0 #FF8C00, stop: 1 white);}
                      '''

        self.ui.tab_a_runhunter_pgb.setStyleSheet(q_pgb_style)
        self.ui.tab_a_runhunter_pgb.hide()
        self.ui.tab_a_msmax_chb.hide()
        self.ui.tab_b_runbatch_pgb.setStyleSheet(q_pgb_style)
        self.ui.tab_b_runbatch_pgb.hide()
        self.ui.tab_c_runlm_pgb.setStyleSheet(q_pgb_style)
        self.ui.tab_c_runlm_pgb.hide()
        self.ui.tab_b_mutlimode_cmb.hide()
        self.ui.tab_b_maxbatch_lb.hide()
        self.ui.tab_b_maxbatch_spb.setValue(1)
        self.ui.tab_b_maxbatch_spb.hide()
        self.ui.tab_b_maxsubcore_spb.setValue(3)
        self.ui.tab_b_maxsubram_spb.setValue(5)
        self.ui.tab_b_maxsubram_spb.setValue(5)
        self.ui.tabframe.removeTab(3)
        self.ui.tab_c_isotopescoremode_lb.hide()
        self.ui.tab_c_isotopescoremode_cmb.hide()
        self.ui.tab_c_scoremode_lb.hide()
        self.ui.tab_c_scoremode_cmb.hide()
        self.ui.tab_c_parallization_cmb.hide()
        self.ui.tab_c_tag_all_fa_chb.hide()

        # define single worker
        self.single_worker = SingleWorker()
        self.single_thread = QtCore.QThread()
        self.single_worker.moveToThread(self.single_thread)
        self.single_worker.workRequested.connect(self.single_thread.start)
        self.single_thread.started.connect(self.single_worker.run_hunter)
        self.single_worker.finished.connect(self.single_worker_on_finish)
        self.single_worker.info_update.connect(self.single_worker_info_update)

        # define batch worker
        self.batch_worker = BatchWorker()
        self.batch_thread = QtCore.QThread()
        self.batch_worker.moveToThread(self.batch_thread)
        self.batch_worker.workRequested.connect(self.batch_thread.start)
        self.batch_thread.started.connect(self.batch_worker.run_hunter)
        self.batch_worker.finished.connect(self.batch_worker_on_finish)
        self.batch_worker.info_update.connect(self.batch_worker_info_update)

        # define lipidmaster worker
        self.lm_worker = LMWorker()
        self.lm_thread = QtCore.QThread()
        self.lm_worker.moveToThread(self.lm_thread)
        self.lm_worker.workRequested.connect(self.lm_thread.start)
        self.lm_thread.started.connect(self.lm_worker.run_lm_generator)
        self.lm_worker.finished.connect(self.lm_worker_on_finish)
        self.lm_worker.info_update.connect(self.lm_worker_info_update)

        # slots for tab a
        QtCore.QObject.connect(self.ui.tab_a_lipidclass_cmb, QtCore.SIGNAL("currentIndexChanged(const QString&)"),
                               self.a_lipid_class_fa_list)
        QtCore.QObject.connect(self.ui.mode_lcms_rb, QtCore.SIGNAL("clicked()"), self.a_set_lc_mode)
        QtCore.QObject.connect(self.ui.mode_static_rb, QtCore.SIGNAL("clicked()"), self.a_set_static_mode)
        QtCore.QObject.connect(self.ui.tab_a_loadfalist_pb, QtCore.SIGNAL("clicked()"), self.a_load_xlsx)
        QtCore.QObject.connect(self.ui.tab_a_loadscorecfg_pb, QtCore.SIGNAL("clicked()"), self.a_loadscore_xlsx)
        # QtCore.QObject.connect(self.ui.tab_a_launchgen_pb, QtCore.SIGNAL("clicked()"), self.a_go_generator)
        QtCore.QObject.connect(self.ui.tab_a_mzml_pb, QtCore.SIGNAL("clicked()"), self.a_load_mzml)
        QtCore.QObject.connect(self.ui.tab_a_saveimgfolder_pb, QtCore.SIGNAL("clicked()"), self.a_save_img2folder)
        QtCore.QObject.connect(self.ui.tab_a_msmax_chb, QtCore.SIGNAL("clicked()"), self.a_max_ms)
        QtCore.QObject.connect(self.ui.tab_a_sumxlsxpath_pb, QtCore.SIGNAL("clicked()"), self.a_save_output)
        # QtCore.QObject.connect(self.ui.tab_a_runhunter_pb, QtCore.SIGNAL("clicked()"), self.a_run_hunter)
        QtCore.QObject.connect(self.ui.tab_a_runhunter_pb, QtCore.SIGNAL("clicked()"), self.single_worker_hunter)
        QtCore.QObject.connect(self.ui.tab_a_cfgpath_pb, QtCore.SIGNAL("clicked()"), self.a_save_cfg)
        QtCore.QObject.connect(self.ui.tab_a_gencfg_pb, QtCore.SIGNAL("clicked()"), self.a_create_cfg)
        # # slots for tab b
        self.ui.tab_b_mutlimode_cmb.currentIndexChanged['QString'].connect(self.b_set_multi_mode)
        QtCore.QObject.connect(self.ui.tab_b_addcfg_pb, QtCore.SIGNAL("clicked()"), self.b_load_batchcfg)
        QtCore.QObject.connect(self.ui.tab_b_addcfgfolder_pb, QtCore.SIGNAL("clicked()"), self.b_load_batchcfgfolder)
        QtCore.QObject.connect(self.ui.tab_b_clearall_pb, QtCore.SIGNAL("clicked()"), self.ui.tab_b_infiles_pte.clear)
        # QtCore.QObject.connect(self.ui.tab_b_runbatch_pb, QtCore.SIGNAL("clicked()"), self.b_run_batchmode)
        QtCore.QObject.connect(self.ui.tab_b_runbatch_pb, QtCore.SIGNAL("clicked()"), self.batch_worker_hunter)
        # # slots for tab c
        QtCore.QObject.connect(self.ui.tab_c_falistpl_pb, QtCore.SIGNAL("clicked()"), self.c_load_falist_pl)
        QtCore.QObject.connect(self.ui.tab_c_falisttg_pb, QtCore.SIGNAL("clicked()"), self.c_load_falist_tg)
        QtCore.QObject.connect(self.ui.tab_c_lmcalcfalist_pb, QtCore.SIGNAL("clicked()"), self.c_load_falist_dg)
        QtCore.QObject.connect(self.ui.tab_c_hgcfg_pb, QtCore.SIGNAL("clicked()"), self.c_load_hgcfg)
        QtCore.QObject.connect(self.ui.tab_c_scorecfgpl_pb, QtCore.SIGNAL("clicked()"), self.c_load_scorecfg_pl)
        QtCore.QObject.connect(self.ui.tab_c_scorecfgtg_pb, QtCore.SIGNAL("clicked()"), self.c_load_scorecfg_tg)
        QtCore.QObject.connect(self.ui.tab_c_scorecfgdg_pb, QtCore.SIGNAL("clicked()"), self.c_load_scorecfg_dg)
        QtCore.QObject.connect(self.ui.tab_c_savesettings_pb, QtCore.SIGNAL("clicked()"), self.c_set_default_cfg)
        QtCore.QObject.connect(self.ui.tab_c_lmexport_pb, QtCore.SIGNAL("clicked()"), self.c_lmexport)
        QtCore.QObject.connect(self.ui.tab_c_lmrun_pb, QtCore.SIGNAL("clicked()"), self.lm_worker_hunter)

        # load configurations

    def load_cfg(self):

        # set click able external links
        for _link_lb in (self.ui.link_gplv2_lb, self.ui.link_source_lb, self.ui.link_tutorial_lb,
                         self.ui.link_paper_lb, self.ui.link_otherprojects_lb, self.ui.link_matplotlib_lb,
                         self.ui.link_numpy_lb, self.ui.link_pandas_lb, self.ui.link_pymzml_lb, self.ui.link_pyside_lb,
                         self.ui.link_bmbf_lb, self.ui.link_emed_lb, self.ui.link_sysmedos_lb, self.ui.link_uni_lb):
            _link_lb.setOpenExternalLinks(True)
        # self.ui.logo_lb.setOpenExternalLinks(True)

        config = configparser.ConfigParser()
        config.read('config.ini')
        if config.has_section('settings'):
            user_cfg = 'settings'
        else:
            if config.has_section('default'):
                user_cfg = 'default'
            else:
                user_cfg = ''
        if len(user_cfg) > 2:
            options = config.options(user_cfg)
            if 'fa_white_list_cfg_pl' in options:
                pl_fawhitelist_path_str = config.get(user_cfg, 'fa_white_list_cfg_pl')
                pl_fawhitelist_path_str, error_log = self.check_file(pl_fawhitelist_path_str, 'FA whitelist for PL')
                if error_log is not None:
                    self.ui.tab_c_falistpl_le.setText(error_log)
                else:
                    self.ui.tab_c_falistpl_le.setText(pl_fawhitelist_path_str)
            if 'lipid_specific_cfg' in options:
                lipid_specific_path_str = config.get(user_cfg, 'lipid_specific_cfg')
                lipid_specific_path_str, error_log = self.check_file(lipid_specific_path_str, 'lipid_specific_cfg')
                if error_log is not None:
                    self.ui.tab_c_hgcfg_le.setText(error_log)
                else:
                    self.ui.tab_c_hgcfg_le.setText(lipid_specific_path_str)
            if 'score_cfg_lpl' in options:
                lpl_score_cfg_path_str = config.get(user_cfg, 'score_cfg_lpl')
                lpl_score_cfg_path_str, error_log = self.check_file(lpl_score_cfg_path_str, 'Score cfg for LPL')
                if error_log is not None:
                    self.ui.tab_c_falisttg_le.setText(error_log)
                else:
                    self.ui.tab_c_falisttg_le.setText(lpl_score_cfg_path_str)
            if 'score_cfg_pl' in options:
                self.ui.tab_c_scorecfgpl_le.setText(config.get(user_cfg, 'score_cfg_pl'))
                pl_score_cfg_path_str = config.get(user_cfg, 'score_cfg_pl')
                pl_score_cfg_path_str, error_log = self.check_file(pl_score_cfg_path_str, 'Score cfg for PL')
                if error_log is not None:
                    self.ui.tab_c_scorecfgpl_le.setText(error_log)
                else:
                    self.ui.tab_c_scorecfgpl_le.setText(pl_score_cfg_path_str)
            if 'score_cfg_tg' in options:
                pl_score_cfg_path_str = config.get(user_cfg, 'score_cfg_tg')
                pl_score_cfg_path_str, error_log = self.check_file(pl_score_cfg_path_str, 'Score cfg for TG')
                if error_log is not None:
                    self.ui.tab_c_scorecfgtg_le.setText(error_log)
                else:
                    self.ui.tab_c_scorecfgtg_le.setText(pl_score_cfg_path_str)
            if 'score_cfg_dg' in options:
                pl_score_cfg_path_str = config.get(user_cfg, 'score_cfg_dg')
                pl_score_cfg_path_str, error_log = self.check_file(pl_score_cfg_path_str, 'Score cfg for DG')
                if error_log is not None:
                    self.ui.tab_c_scorecfgdg_le.setText(error_log)
                else:
                    self.ui.tab_c_scorecfgdg_le.setText(pl_score_cfg_path_str)
            if 'score_mode' in options:
                if config.get(user_cfg, 'score_mode').upper() in ['RANK', '']:
                    self.ui.tab_c_scoremode_cmb.setCurrentIndex(0)
                elif config.get(user_cfg, 'score_mode') is None:
                    self.ui.tab_c_scoremode_cmb.setCurrentIndex(0)
                elif config.get(user_cfg, 'score_mode').upper() == 'INTENSITY':
                    self.ui.tab_c_scoremode_cmb.setCurrentIndex(1)
                else:
                    self.ui.tab_c_scoremode_cmb.setCurrentIndex(0)
            else:
                self.ui.tab_c_scoremode_cmb.setCurrentIndex(0)
            if 'isotope_13c_mode' in options:
                if config.get(user_cfg, 'isotope_13c_mode').upper() == 'ON':
                    self.ui.tab_c_isotopescoremode_cmb.setCurrentIndex(1)
                elif config.get(user_cfg, 'isotope_13c_mode').upper() == 'OFF':
                    self.ui.tab_c_isotopescoremode_cmb.setCurrentIndex(0)
                else:
                    self.ui.tab_c_isotopescoremode_cmb.setCurrentIndex(0)
            else:
                self.ui.tab_c_isotopescoremode_cmb.setCurrentIndex(0)
            if 'parallel_target' in options:
                if config.get(user_cfg, 'parallel_target').upper() == 'CPU':
                    self.ui.tab_c_parallization_cmb.setCurrentIndex(0)
                elif config.get(user_cfg, 'parallel_target') in ['CPU_and_GPU', 'GPU', 'CPUandGPU', 'CPUGPU',
                                                                 'parallel']:
                    self.ui.tab_c_parallization_cmb.setCurrentIndex(1)
                else:
                    self.ui.tab_c_parallization_cmb.setCurrentIndex(0)
            if 'max_cpu_core' in options:
                self.ui.tab_c_cores_spb.setValue(int(config.get(user_cfg, 'max_cpu_core')))
            if 'max_ram' in options:
                self.ui.tab_c_ram_spb.setValue(int(config.get(user_cfg, 'max_ram')))
            if 'img_type' in options:
                if config.get(user_cfg, 'img_type').upper() == 'PNG':
                    self.ui.tab_c_imagetype_cmb.setCurrentIndex(0)
                elif config.get(user_cfg, 'img_type').upper() == 'SVG':
                    self.ui.tab_c_imagetype_cmb.setCurrentIndex(1)
                else:
                    self.ui.tab_c_imagetype_cmb.setCurrentIndex(0)
            if 'img_dpi' in options:
                self.ui.tab_c_dpi_spb.setValue(int(config.get(user_cfg, 'img_dpi')))
            else:
                self.ui.tab_c_dpi_spb.setValue(300)
            if 'tag_all_fa_check' in options:
                if config.get(user_cfg, 'tag_all_fa_check').upper() == 'ON':
                    self.ui.tab_c_tag_all_fa_chb.setChecked(True)
                else:
                    self.ui.tab_c_tag_all_fa_chb.setChecked(False)
            else:
                self.ui.tab_c_tag_all_fa_chb.setChecked(True)
            # UI settings
            if 'main_tab' in options and 1 <= int(config.get(user_cfg, 'main_tab')) <= 5:
                self.ui.tabframe.setCurrentIndex(int(config.get(user_cfg, 'main_tab')) - 1)
            else:
                self.ui.tabframe.setCurrentIndex(0)
            if 'runhunter_tab' in options and 1 <= int(config.get(user_cfg, 'runhunter_tab')) <= 2:
                self.ui.runhunter_tabframe.setCurrentIndex(int(config.get(user_cfg, 'runhunter_tab')) - 1)
            else:
                self.ui.runhunter_tabframe.setCurrentIndex(0)
            if 'lipidgen_tab' in options and 1 <= int(config.get(user_cfg, 'lipidgen_tab')) <= 2:
                self.ui.lipidgen_tabframe.setCurrentIndex(int(config.get(user_cfg, 'lipidgen_tab')) - 1)
            else:
                self.ui.lipidgen_tabframe.setCurrentIndex(0)

            if 'default_vendor' in options and 0 <= int(config.get(user_cfg, 'default_vendor')) <= 4:
                self.ui.vendor_cmb.setCurrentIndex(int(config.get(user_cfg, 'default_vendor')))
            else:
                self.ui.vendor_cmb.setCurrentIndex(0)
            if 'default_lipid' in options and 0 <= int(config.get(user_cfg, 'default_lipid')) <= 9:
                self.ui.tab_a_lipidclass_cmb.setCurrentIndex(int(config.get(user_cfg, 'default_lipid')))
            else:
                self.ui.tab_a_lipidclass_cmb.setCurrentIndex(0)

    @staticmethod
    def get_same_files(folder, filetype_lst):
        """
        find all files with same type in specified folder
        :param str folder: absolute file path
        :param list filetype_lst: e.g. ['*.mzml', '*.mzML']
        :return: a list of absolute file path
        :rtype: list
        """
        if folder is not '':
            os.chdir(folder)
            _pre_found_lst = []
            for _filetype in filetype_lst:
                _tmp_found_lst = glob.glob(_filetype)
                # merge list
                _pre_found_lst += [f for f in _tmp_found_lst if f not in _pre_found_lst]
            filename_lst = _pre_found_lst
            abs_path_lst = list(os.path.abspath(ff) for ff in _pre_found_lst)
        else:
            filename_lst = []
            abs_path_lst = []

        return filename_lst, abs_path_lst

    def open_file(self, info_str, lb_obj):
        open_file_dialog = QtGui.QFileDialog(self)
        open_file_dialog.setNameFilters([info_str])
        open_file_dialog.selectNameFilter(info_str)
        if open_file_dialog.exec_():
            lb_obj.clear()
            file_str = open_file_dialog.selectedFiles()[0]
            file_str = os.path.abspath(file_str)
            lb_obj.setText(file_str)

    @staticmethod
    def check_file(usr_path, info_str):
        file_abs_path = ''
        try:
            if os.path.isfile(usr_path):
                error_log = None
                file_abs_path = os.path.abspath(usr_path)
            else:
                error_log = '!! Failed to load {_file} !!'.format(_file=info_str)
        except IOError:
            error_log = '!! Failed to load {_file} !!'.format(_file=info_str)

        return file_abs_path, error_log

    @staticmethod
    def check_folder(usr_path, info_str):
        folder_abs_path = ''
        try:
            if os.path.isdir(usr_path):
                print('Folder existed...\n', usr_path)
                error_log = None
                folder_abs_path = os.path.abspath(usr_path)
                print('abs path of folder\n', folder_abs_path)
            # else:
            #     if platform == "linux" or platform == "linux2":
            #         l_cwd = os.getcwd()
            #         os.chdir('/')
            #         if os.path.isdir(usr_path):
            #             print('Folder existed...\n', usr_path)
            #             error_log = None
            #             folder_abs_path = os.path.abspath(usr_path)
            #             print('abs path of folder\n', folder_abs_path)
            #         else:
            #             if os.path.isdir('/' + usr_path):
            #                 print('Folder existed...\n', usr_path)
            #                 error_log = None
            #                 folder_abs_path = os.path.abspath(usr_path)
            #                 print('abs path of folder\n', folder_abs_path)
            #             else:
            #                 print('No folder...\n', usr_path)
            #                 os.makedirs(usr_path)
            #                 print('Folder created... %s' % usr_path)
            #                 error_log = ''
            #         os.chdir(l_cwd)
            else:
                print('No folder...\n', usr_path)
                os.makedirs(usr_path)
                print('Folder created... %s' % usr_path)
                error_log = ''
        except IOError:
            error_log = '!! Failed to open folder {_file} !!'.format(_file=info_str)

        return folder_abs_path, error_log

    def a_load_xlsx(self):
        file_info_str = 'FA white list files (*.xlsx *.XLSX)'
        self.open_file(file_info_str, self.ui.tab_a_loadfalist_le)

    def a_loadscore_xlsx(self):
        file_info_str = 'Weight list files (*.xlsx *.XLSX)'
        self.open_file(file_info_str, self.ui.tab_a_loadscorecfg_le)

    # def a_go_generator(self):
    #     self.ui.tabframe.setCurrentIndex(3)

    def a_load_mzml(self):
        file_info_str = 'mzML spectra files (*.mzML *.mzml)'
        self.open_file(file_info_str, self.ui.tab_a_mzml_le)

    def a_save_img2folder(self):
        a_save_img2folder_str = QtGui.QFileDialog.getExistingDirectory()
        self.ui.tab_a_saveimgfolder_le.clear()
        self.ui.tab_a_saveimgfolder_le.setText(a_save_img2folder_str)

    def a_save_output(self):
        a_save_output_path = QtGui.QFileDialog.getSaveFileName(caption='Save file', filter='.xlsx')
        self.ui.tab_a_savexlsxpath_le.clear()
        usr_output_path = a_save_output_path[0]
        if usr_output_path[-5:] != '.xlsx':
            usr_output_path += '.xlsx'
        a_save_output_str = os.path.abspath(usr_output_path)

        self.ui.tab_a_savexlsxpath_le.setText(a_save_output_str)

    def a_max_ms(self):
        if self.ui.tab_a_msmax_chb.isChecked():
            self.ui.tab_a_msmax_spb.show()
            self.ui.tab_a_msmax_spb.setValue(100 * self.ui.tab_a_msthreshold_spb.value())
        else:
            self.ui.tab_a_msmax_spb.hide()
            self.ui.tab_a_msmax_spb.setValue(0)

    def a_set_lc_mode(self):
        self.ui.tab_a_rtstart_dspb.setValue(10.0)
        self.ui.tab_a_rtend_dspb.setValue(25.0)

    def a_set_static_mode(self):
        self.ui.tab_a_rtstart_dspb.setValue(0.1)
        self.ui.tab_a_rtend_dspb.setValue(15.0)

    def a_lipid_class_fa_list(self):
        _lipid_class_info = str(self.ui.tab_a_lipidclass_cmb.currentText())
        lipid_class_checker = re.compile(r'(.*)( [(])(\w{2,3})([)] )(.*)')
        lipid_class_match = lipid_class_checker.match(_lipid_class_info)

        if lipid_class_match:
            lipid_class_info_lst = lipid_class_match.groups()
            _lipid_class = lipid_class_info_lst[2]
            _lipid_charge = lipid_class_info_lst[4]
        else:
            _lipid_class = ''
            _lipid_charge = ''

        pl_fa_cfg = self.ui.tab_c_falistpl_le.text()
        lpl_score_cfg = self.ui.tab_c_falisttg_le.text()
        pl_score_cfg = self.ui.tab_c_scorecfgpl_le.text()
        tg_score_cfg = self.ui.tab_c_scorecfgtg_le.text()
        dg_score_cfg = self.ui.tab_c_scorecfgdg_le.text()
        usr_fa_cfg = self.ui.tab_a_loadfalist_le.text()
        usr_score_cfg = self.ui.tab_a_loadscorecfg_le.text()

        if _lipid_class in ['PA', 'PC', 'PE', 'PG', 'PI', 'PIP', 'PS']:
            # if usr_fa_cfg in ['', tg_fa_cfg]:
            self.ui.tab_a_loadfalist_le.setText(pl_fa_cfg)
            # if usr_score_cfg in ['', dg_score_cfg, tg_score_cfg]:
            self.ui.tab_a_loadscorecfg_le.setText(pl_score_cfg)
            self.ui.tab_a_mzstart_dspb.setValue(600)
            self.ui.tab_a_mzend_dspb.setValue(1000)
            self.ui.tab_a_score_spb.setValue(40)
            if self.ui.mode_lcms_rb.isChecked():
                self.ui.tab_a_rtstart_dspb.setValue(10.0)
                self.ui.tab_a_rtend_dspb.setValue(25.0)
        elif _lipid_class in ['LPA', 'LPC', 'LPE', 'LPG', 'LPI', 'LPIP', 'LPS']:
            # if usr_fa_cfg in ['', tg_fa_cfg]:
            self.ui.tab_a_loadfalist_le.setText(pl_fa_cfg)
            # if usr_score_cfg in ['', dg_score_cfg, tg_score_cfg]:
            self.ui.tab_a_loadscorecfg_le.setText(lpl_score_cfg)
            self.ui.tab_a_mzstart_dspb.setValue(300)
            self.ui.tab_a_mzend_dspb.setValue(800)
            self.ui.tab_a_score_spb.setValue(25)
            if self.ui.mode_lcms_rb.isChecked():
                self.ui.tab_a_rtstart_dspb.setValue(3.0)
                self.ui.tab_a_rtend_dspb.setValue(15.0)
        elif _lipid_class in ['TG', 'DG', 'MG']:
            # if usr_fa_cfg in ['', pl_fa_cfg]:
            self.ui.tab_a_loadfalist_le.setText(pl_fa_cfg)
            self.ui.tab_a_mzstart_dspb.setValue(600)
            self.ui.tab_a_mzend_dspb.setValue(1200)
            self.ui.tab_a_score_spb.setValue(50)
            if self.ui.mode_lcms_rb.isChecked():
                self.ui.tab_a_rtstart_dspb.setValue(20.0)
                self.ui.tab_a_rtend_dspb.setValue(30.0)
            if _lipid_class in ['TG'] and _lipid_charge not in ['[M+Na]+']:
                # if usr_score_cfg in ['', pl_score_cfg, dg_score_cfg]:
                self.ui.tab_a_loadscorecfg_le.setText(tg_score_cfg)
            elif _lipid_class in ['TG'] and _lipid_charge in ['[M+Na]+']:
                tg_na_abs_path, error_log = self.check_file('.\ConfigurationFiles\\2-Score_weight_TG_Na.xlsx',
                                                            'TG [M+Na]+ Weight factor'
                                                            )
                self.ui.tab_a_loadscorecfg_le.setText(tg_na_abs_path)
            elif _lipid_class in ['DG']:
                # if usr_score_cfg in ['', pl_score_cfg, tg_score_cfg]:
                self.ui.tab_a_loadscorecfg_le.setText(dg_score_cfg)
                self.ui.tab_a_mzstart_dspb.setValue(300)
                self.ui.tab_a_mzend_dspb.setValue(900)
                self.ui.tab_a_score_spb.setValue(40)
                if self.ui.mode_lcms_rb.isChecked():
                    self.ui.tab_a_rtstart_dspb.setValue(10.0)
                    self.ui.tab_a_rtend_dspb.setValue(20.0)
            else:
                self.ui.tab_a_loadscorecfg_le.setText('')
        else:
            self.ui.tab_a_loadfalist_le.setText('')
            self.ui.tab_a_loadscorecfg_le.setText('')

    def a_get_params(self):

        error_log_lst = []

        usr_vendor_str = self.ui.vendor_cmb.currentText().lower()
        vendors_dct = {'therm': 'thermo', 'water': 'waters', 'sciex': 'sciex', 'agile': 'agilent'}
        if usr_vendor_str[0:5] in list(vendors_dct.keys()):
            usr_vendor = vendors_dct[usr_vendor_str[0:5]]
        else:
            usr_vendor = ''
            error_log_lst.append('!! Please select an instrument vendor!!')

        if self.ui.mode_lcms_rb.isChecked():
            usr_exp_mode = 'LC-MS'
        elif self.ui.mode_static_rb.isChecked():
            usr_exp_mode = 'Shotgun'
        else:
            usr_exp_mode = 'LC-MS'

        _lipid_class_info = str(self.ui.tab_a_lipidclass_cmb.currentText())

        lipid_class_checker = re.compile(r'(.*)( [(])(\w{2,3})([)] )(.*)')

        lipid_class_match = lipid_class_checker.match(_lipid_class_info)

        if lipid_class_match:
            lipid_class_info_lst = lipid_class_match.groups()
            _lipid_class = lipid_class_info_lst[2]
            _lipid_charge = lipid_class_info_lst[4]
        else:
            _lipid_class = ''
            _lipid_charge = ''
            error_log_lst.append('!! Please select a lipid class!!')

        # if _lipid_class in ['PA', 'PC', 'PE', 'PG', 'PI', 'PIP', 'PS']:
        #     score_cfg = self.ui.tab_a_loadscorecfg_le.text()
        # elif _lipid_class in ['TG', 'DG', 'MG']:
        #     score_cfg = self.ui.tab_a_loadscorecfg_le.text()
        # else:
        #     score_cfg = ''
        #     error_log_lst.append('!! No Corresponding Score weight factor for selected lipid class!!')

        score_cfg = self.ui.tab_a_loadscorecfg_le.text()
        fawhitelist_path_str = str(self.ui.tab_a_loadfalist_le.text())
        mzml_path_str = str(self.ui.tab_a_mzml_le.text())
        img_output_folder_str = str(self.ui.tab_a_saveimgfolder_le.text())
        if img_output_folder_str[-1] in ['\\', '/']:
            img_output_folder_str = img_output_folder_str[:-1]
            self.ui.tab_a_saveimgfolder_le.clear()
            self.ui.tab_a_saveimgfolder_le.setText(img_output_folder_str)
        else:
            pass
        xlsx_output_path_str = str(self.ui.tab_a_savexlsxpath_le.text())
        if xlsx_output_path_str[-5:] != '.xlsx':
            xlsx_output_path_str += '.xlsx'
            self.ui.tab_a_savexlsxpath_le.clear()
            self.ui.tab_a_savexlsxpath_le.setText(xlsx_output_path_str)

        rt_start = self.ui.tab_a_rtstart_dspb.value()
        rt_end = self.ui.tab_a_rtend_dspb.value()
        mz_start = self.ui.tab_a_mzstart_dspb.value()
        mz_end = self.ui.tab_a_mzend_dspb.value()
        dda_top = self.ui.tab_a_dda_spb.value()
        pr_window = self.ui.tab_a_prwindow_spb.value()
        isotope_score_filter = self.ui.tab_a_isotopescore_spb.value()
        rank_score_filter = self.ui.tab_a_score_spb.value()
        ms_th = self.ui.tab_a_msthreshold_spb.value()
        ms2_th = self.ui.tab_a_ms2threshold_spb.value()
        ms_ppm = self.ui.tab_a_msppm_spb.value()
        ms2_ppm = self.ui.tab_a_ms2ppm_spb.value()
        ms2_info_threshold = self.ui.tab_a_ms2infoth_dspb.value() * 0.01

        ms_max = 0
        if self.ui.tab_a_msmax_chb.isChecked() and self.ui.tab_a_msmax_spb.value() > ms_th + 1:
            ms_max = self.ui.tab_a_msmax_spb.value()
        lipid_specific_cfg = self.ui.tab_c_hgcfg_le.text()

        core_num = self.ui.tab_c_cores_spb.value()
        max_ram = self.ui.tab_c_ram_spb.value()
        img_typ = self.ui.tab_c_imagetype_cmb.currentText()[1:]
        img_dpi = self.ui.tab_c_dpi_spb.value()

        fawhitelist_path_str, error_log = self.check_file(fawhitelist_path_str, 'FA whitelist')
        if error_log is not None:
            error_log_lst.append(error_log)
        lipid_specific_cfg, error_log = self.check_file(lipid_specific_cfg, 'configuration for Phospholipids')
        if error_log is not None:
            error_log_lst.append(error_log)
        score_cfg, error_log = self.check_file(score_cfg, 'W_frag score configuration')
        if error_log is not None:
            error_log_lst.append(error_log)
        abs_img_output_folder_str, error_log = self.check_folder(img_output_folder_str, 'Output folder')
        if error_log is not None:
            error_log_lst.append(error_log)
            self.ui.tab_a_saveimgfolder_le.setText(error_log)
        else:
            # if platform == "linux" or platform == "linux2":
            #     pass
            # else:
            if abs_img_output_folder_str != img_output_folder_str:
                self.ui.tab_a_saveimgfolder_le.clear()
                self.ui.tab_a_saveimgfolder_le.setText(abs_img_output_folder_str)
                print('!! Image output folder not correct !!')
                print('>> Propose to save in folder: %s' % abs_img_output_folder_str)
                img_output_folder_str = abs_img_output_folder_str

        if xlsx_output_path_str[-5:] == '.xlsx':
            pass
        else:
            xlsx_output_path_str = '%s.xlsx' % xlsx_output_path_str
            self.ui.tab_a_savexlsxpath_le.setText(xlsx_output_path_str)
            error_log_lst.append('!! Excel output path not correct !!')
            error_log_lst.append('>> Propose to save as: %s' % xlsx_output_path_str)

        usr_score_mode = self.ui.tab_c_scoremode_cmb.currentIndex()
        if usr_score_mode == 0:
            rank_score = True
        else:
            rank_score = False

        usr_isotope_score_mode = self.ui.tab_c_isotopescoremode_cmb.currentIndex()
        if usr_isotope_score_mode == 0:
            fast_isotope = False
        else:
            fast_isotope = True

        if self.ui.tab_c_tag_all_fa_chb.isChecked():
            tag_all_sn = True
        else:
            tag_all_sn = False

        start_time_str = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())

        score_filter = rank_score_filter

        hunter_param_dct = {'vendor': usr_vendor, 'experiment_mode': usr_exp_mode,
                            'lipid_class': _lipid_class, 'charge_mode': _lipid_charge,
                            'fawhitelist_path_str': fawhitelist_path_str,
                            'score_cfg': score_cfg,
                            'mzml_path_str': mzml_path_str,
                            'img_output_folder_str': img_output_folder_str,
                            'xlsx_output_path_str': xlsx_output_path_str,
                            'rt_start': rt_start, 'rt_end': rt_end,
                            'mz_start': mz_start, 'mz_end': mz_end,
                            'dda_top': dda_top, 'pr_window': pr_window,
                            'ms_th': ms_th, 'ms_ppm': ms_ppm,
                            'ms2_th': ms2_th, 'ms2_ppm': ms2_ppm, 'ms2_infopeak_threshold': ms2_info_threshold,
                            'rank_score_filter': rank_score_filter, 'score_filter': score_filter,
                            'isotope_score_filter': isotope_score_filter,
                            'lipid_specific_cfg': lipid_specific_cfg,
                            'core_number': core_num, 'max_ram': max_ram,
                            'img_type': img_typ, 'img_dpi': img_dpi,
                            'hunter_folder': self.lipidhunter_cwd,
                            'hunter_start_time': start_time_str,
                            }
        debug_mode = 'OFF'
        if debug_mode == 'ON':
            hunter_param_dct['rank_score'] = rank_score
            hunter_param_dct['tag_all_sn'] = tag_all_sn
            hunter_param_dct['fast_isotope'] = fast_isotope
            hunter_param_dct['ms_max'] = ms_max
        else:
            hunter_param_dct['rank_score'] = True
            hunter_param_dct['tag_all_sn'] = True
            hunter_param_dct['fast_isotope'] = False
            hunter_param_dct['ms_max'] = 0

            return hunter_param_dct, error_log_lst

    def a_save_cfg(self):
        a_save_cfg_path = QtGui.QFileDialog.getSaveFileName(caption='Save file', filter='.txt')
        self.ui.tab_a_cfgpath_le.clear()
        usr_cfg_path = a_save_cfg_path[0]
        if usr_cfg_path[-4:] != '.txt':
            usr_cfg_path += '.txt'
        a_save_cfg_str = os.path.abspath(usr_cfg_path)
        self.ui.tab_a_cfgpath_le.setText(a_save_cfg_str)

    def a_create_cfg(self):

        param_cfg_path_str = str(self.ui.tab_a_cfgpath_le.text())
        try:
            if param_cfg_path_str[-4:] == '.txt':
                pass
            else:
                param_cfg_path_str = '%s.txt' % param_cfg_path_str
                self.ui.tab_a_cfgpath_le.clear()
                self.ui.tab_a_cfgpath_le.setText(param_cfg_path_str)
            param_cfg_directory = os.path.dirname(param_cfg_path_str)
            if not os.path.exists(param_cfg_directory):
                os.makedirs(param_cfg_directory)
            abs_param_cfg_directory_str, error_log = self.check_folder(param_cfg_directory, 'Settings for batch mode')
            if error_log is not None:
                self.ui.tab_a_gencfg_pte.insertPlainText(error_log)
            if abs_param_cfg_directory_str != param_cfg_directory:
                param_cfg_path = os.path.split(param_cfg_path_str)
                param_cfg_path_str = os.path.join(abs_param_cfg_directory_str, param_cfg_path[1])
                self.ui.tab_a_gencfg_pte.insertPlainText('>>> try to save as: %s\n' % param_cfg_path_str)

        except Exception as _err:
            print(_err)
            self.ui.tab_a_gencfg_pte.insertPlainText('!!! Failed to save settings for batch mode !!!\n')
            self.ui.tab_a_gencfg_pte.insertPlainText('!! Can not save settings to file !!\n')
            self.ui.tab_a_gencfg_pte.insertPlainText(str(_err))

        hunter_param_dct, error_log_lst = self.a_get_params()
        error_log_lst = [_f for _f in error_log_lst if _f]

        if len(error_log_lst) > 0:
            print('Parameter error:', error_log_lst)
            error_log_lst.append('!!! Failed to save settings for batch mode !!!')
            error_log_lst.append('>>> Please check your settings, and try to save again ...\n')
            self.ui.tab_a_gencfg_pte.appendPlainText('\n'.join(error_log_lst) + '\n')
        else:

            try:

                config = configparser.ConfigParser()
                with open(param_cfg_path_str, 'w') as usr_param_cfg:
                    config.add_section('parameters')
                    for param in list(hunter_param_dct.keys()):
                        config.set('parameters', str(param), str(hunter_param_dct[param]))
                    config.write(usr_param_cfg)
                    self.ui.tab_a_gencfg_pte.insertPlainText('>>> Successfully saved as: \n')
                    self.ui.tab_a_gencfg_pte.insertPlainText(param_cfg_path_str)
                    self.ui.tab_a_gencfg_pte.insertPlainText('\n')
            except Exception as _err:
                print(_err)
                self.ui.tab_a_gencfg_pte.insertPlainText(str(_err))

    def b_set_multi_mode(self):

        multi_mode_idx = self.ui.tab_b_mutlimode_cmb.currentIndex()

        if multi_mode_idx == 1:
            print('Set Batch mode to: Multi processing mode')
            self.ui.tab_b_maxbatch_lb.show()
            self.ui.tab_b_maxbatch_spb.show()
            self.ui.tab_b_maxsubcore_lb.show()
            self.ui.tab_b_maxsubcore_spb.show()
            self.ui.tab_b_maxsubram_lb.show()
            self.ui.tab_b_maxsubram_spb.show()
        elif multi_mode_idx == 0:
            print('Set Batch mode to: Single processing mode')
            self.ui.tab_b_maxbatch_lb.hide()
            self.ui.tab_b_maxbatch_spb.hide()
            self.ui.tab_b_maxsubcore_lb.hide()
            self.ui.tab_b_maxsubcore_spb.hide()
            self.ui.tab_b_maxsubram_lb.hide()
            self.ui.tab_b_maxsubram_spb.hide()

    @staticmethod
    def b_get_same_files(folder, filetype_lst):

        if folder is not '':
            os.chdir(folder)
            _pre_found_lst = []
            for _filetype in filetype_lst:
                _tmp_found_lst = glob.glob(_filetype)
                # merge list
                _pre_found_lst += [f for f in _tmp_found_lst if f not in _pre_found_lst]
            filename_lst = _pre_found_lst
            abs_path_lst = list(os.path.abspath(ff) for ff in _pre_found_lst)
        else:
            filename_lst = []
            abs_path_lst = []

        return filename_lst, abs_path_lst

    def b_load_batchcfg(self):
        # check existed files
        _loaded_files = str(self.ui.tab_b_infiles_pte.toPlainText())
        _loaded_lst = _loaded_files.split('\n')

        b_load_cfg_dialog = QtGui.QFileDialog(self)
        b_load_cfg_dialog.setNameFilters(['LipidHunter batch mode files (*.txt)'])
        b_load_cfg_dialog.selectNameFilter('LipidHunter batch mode files (*.txt)')
        if b_load_cfg_dialog.exec_():
            b_load_cfg_str = b_load_cfg_dialog.selectedFiles()[0]
            b_load_cfg_str = os.path.abspath(b_load_cfg_str)
            if b_load_cfg_str not in _loaded_lst:
                self.ui.tab_b_infiles_pte.insertPlainText(b_load_cfg_str)  # take unicode only
                self.ui.tab_b_infiles_pte.insertPlainText('\n')
            else:
                _msgBox = QtGui.QMessageBox()
                _msgBox.setText('Batch config file has been chosen already.')
                _msgBox.exec_()

    def b_load_batchcfgfolder(self):
        # check existed files
        _loaded_files = str(self.ui.tab_b_infiles_pte.toPlainText())
        _loaded_lst = _loaded_files.split('\n')

        b_load_cfgfolder_str = QtGui.QFileDialog.getExistingDirectory()
        _cfg_name_lst, _cfg_path_lst = self.b_get_same_files(b_load_cfgfolder_str, filetype_lst=['*.txt', '*.txt'])
        _duplicated_str = ''
        for _cfg in _cfg_path_lst:
            if _cfg not in _loaded_lst:
                self.ui.tab_b_infiles_pte.insertPlainText(_cfg)
                self.ui.tab_b_infiles_pte.insertPlainText('\n')
            else:
                _duplicated_str = _duplicated_str + _cfg + '\n'
        if len(_duplicated_str) > 0:
            _msgBox = QtGui.QMessageBox()
            _msgBox.setText(_duplicated_str + 'Already chosen. \n Skipped')
            _msgBox.exec_()

    @staticmethod
    def b_read_cfg(batch_cfg):
        cfg_params_dct = {}
        cfg_error = ''

        i_type_key_lst = ['ms_th', 'ms2_th', 'ms_ppm', 'ms2_ppm', 'dda_top', 'sn_ratio',
                          'core_number', 'max_ram', 'img_dpi', 'ms_max']
        f_type_key_lst = ['rt_start', 'rt_end', 'mz_start', 'mz_end', 'pr_window', 'ms2_infopeak_threshold',
                          'score_filter', 'isotope_score_filter', 'rank_score_filter']
        b_type_key_lst = ['rank_score', 'fast_isotope', 'tag_all_sn']

        print('Input LipidHunter configuration file : ', batch_cfg)
        if os.path.isfile(batch_cfg):
            with open(batch_cfg) as _cfg_obj:
                config = configparser.ConfigParser()
                try:
                    try:
                        config.read_file(_cfg_obj)
                    except AttributeError:  # for python 2.7.14
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
                except Exception as _cfg_e:
                    print('ERROR', _cfg_e)
                    cfg_error += str(_cfg_e)
                    cfg_error += '\n'
        return cfg_params_dct, cfg_error

    def c_set_default_cfg(self):
        config = configparser.ConfigParser()
        with open('config.ini', 'w') as default_cfg:
            config.add_section('settings')
            config.set('settings', 'fa_white_list_cfg_pl', self.ui.tab_c_falistpl_le.text())
            config.set('settings', 'fa_white_list_cfg_tg', self.ui.tab_c_falisttg_le.text())
            config.set('settings', 'score_cfg_pl', self.ui.tab_c_scorecfgpl_le.text())
            config.set('settings', 'score_cfg_tg', self.ui.tab_c_scorecfgtg_le.text())
            config.set('settings', 'score_cfg_dg', self.ui.tab_c_scorecfgdg_le.text())
            config.set('settings', 'lipid_specific_cfg', self.ui.tab_c_hgcfg_le.text())

            if self.ui.tab_c_scoremode_cmb.currentIndex() == 0:
                config.set('settings', 'score_mode', 'RANK')
            elif self.ui.tab_c_scoremode_cmb.currentIndex() == 1:
                config.set('settings', 'score_mode', 'INTENSITY')
            else:
                config.set('settings', 'score_mode', 'RANK')

            if self.ui.tab_c_isotopescoremode_cmb.currentIndex() == 1:
                config.set('settings', 'isotope_13c_mode', 'ON')
            else:
                config.set('settings', 'isotope_13c_mode', 'OFF')

            # parallel processing settings
            if self.ui.tab_c_parallization_cmb.currentIndex() == 0:
                config.set('settings', 'parallel_target', 'CPU')
            elif self.ui.tab_c_parallization_cmb.currentIndex() == 1:
                config.set('settings', 'parallel_target', 'CPU_and_GPU')
            else:
                config.set('settings', 'parallel_target', 'CPU')
            config.set('settings', 'max_cpu_core', str(self.ui.tab_c_cores_spb.value()))
            config.set('settings', 'max_ram', str(self.ui.tab_c_ram_spb.value()))

            # other settings
            config.set('settings', 'img_type', str(self.ui.tab_c_imagetype_cmb.currentText()).upper()[1:])
            config.set('settings', 'img_dpi', str(self.ui.tab_c_dpi_spb.value()))

            if self.ui.tab_c_tag_all_fa_chb.isChecked():
                config.set('settings', 'tag_all_fa_check', 'ON')
            else:
                config.set('settings', 'tag_all_fa_check', 'OFF')

            config.set('settings', 'main_tab', '1')
            config.set('settings', 'runhunter_tab', '1')
            config.set('settings', 'lipidgen_tab', '1')
            config.set('settings', 'default_vendor', '0')
            config.set('settings', 'default_lipid', '0')

            config.write(default_cfg)

    def c_load_falist_pl(self):
        file_info_str = 'FA white list files (*.xlsx *.XLSX)'
        self.open_file(file_info_str, self.ui.tab_c_falistpl_le)

    def c_load_falist_tg(self):
        file_info_str = 'FA white list files (*.xlsx *.XLSX)'
        self.open_file(file_info_str, self.ui.tab_c_falisttg_le)

    def c_load_falist_dg(self):
        file_info_str = 'FA white list files (*.xlsx *.XLSX)'
        self.open_file(file_info_str, self.ui.tab_c_lmcalcfalist_le)

    def c_load_hgcfg(self):
        file_info_str = 'MS Excel files (*.xlsx *.XLSX)'
        self.open_file(file_info_str, self.ui.tab_c_hgcfg_le)

    def c_load_scorecfg_pl(self):
        file_info_str = 'MS Excel files (*.xlsx *.XLSX)'
        self.open_file(file_info_str, self.ui.tab_c_scorecfgpl_le)

    def c_load_scorecfg_tg(self):
        file_info_str = 'MS Excel files (*.xlsx *.XLSX)'
        self.open_file(file_info_str, self.ui.tab_c_scorecfgtg_le)

    def c_load_scorecfg_dg(self):
        file_info_str = 'MS Excel files (*.xlsx *.XLSX)'
        self.open_file(file_info_str, self.ui.tab_c_scorecfgdg_le)

    def c_lmexport(self):
        c_lmexport_path = QtGui.QFileDialog.getSaveFileName(caption='Save file', filter='.csv')
        self.ui.tab_c_lmexport_le.clear()
        usr_lm_path = c_lmexport_path[0]
        if usr_lm_path[-4:] != '.csv':
            usr_lm_path += '.csv'
        c_lmexport_str = os.path.abspath(usr_lm_path)
        self.ui.tab_c_lmexport_le.setText(c_lmexport_str)

    def single_worker_on_finish(self):
        self.single_thread.quit()
        print('!! single_worker stopped !!')
        self.ui.tab_a_runhunter_pb.setText(QtGui.QApplication.translate('MainWindow', 'Hunt for lipids!', None,
                                                                        QtGui.QApplication.UnicodeUTF8))

        self.ui.tab_a_runhunter_pgb.setMinimum(0)
        self.ui.tab_a_runhunter_pgb.setMaximum(100)
        self.ui.tab_a_runhunter_pgb.hide()
        self.ui.tab_a_runhunter_pb.setEnabled(True)
        self.ui.tab_b_runbatch_pb.setEnabled(True)
        self.ui.tab_c_lmrun_pb.setEnabled(True)

    def single_worker_hunter(self):

        self.ui.tab_a_statusrun_pte.clear()
        self.ui.tab_a_statusrun_pte.setPlainText('')

        ready_to_run = False

        hunter_param_dct, error_log_lst = self.a_get_params()

        print('Vendor mode = %s, Experiment mode = %s' % (hunter_param_dct['vendor'],
                                                          hunter_param_dct['experiment_mode']))
        print('Isotope score mode = %s' % (hunter_param_dct['fast_isotope']))
        print('Rankscore mode = %s' % (hunter_param_dct['rank_score']))
        print('Hunter started!')

        output_folder_path = hunter_param_dct['img_output_folder_str']

        if platform == "linux" or platform == "linux2":
            l_cwd = os.getcwd()
            os.chdir('/')
            if os.path.isdir(output_folder_path):
                print('Folder existed...\n', output_folder_path)
                error_log = None
                folder_abs_path = os.path.abspath(output_folder_path)
                print('abs path of folder\n', folder_abs_path)
            else:
                if os.path.isdir('/' + output_folder_path):
                    print('Folder existed...\n', output_folder_path)
                    error_log = None
                    folder_abs_path = os.path.abspath(output_folder_path)
                    print('abs path of folder\n', folder_abs_path)
                else:
                    print('No folder...\n', output_folder_path)
                    os.makedirs(output_folder_path)
                    print('Folder created... %s' % output_folder_path)
                    error_log = ''
            os.chdir(l_cwd)
        else:
            if os.path.isdir(output_folder_path):
                print('Output folder path... %s' % output_folder_path)
            else:
                try:
                    os.mkdir(output_folder_path)
                    print('Output folder created... %s' % output_folder_path)
                except IOError:
                    error_log_lst.append('!! Failed to create output folder !!')

        param_log_output_path_str = os.path.join(output_folder_path + '/LipidHunter_Params-Log_%s.txt'
                                                 % hunter_param_dct['hunter_start_time'])
        print('param_log_output_path_str', param_log_output_path_str)

        try:
            config = configparser.ConfigParser()
            with open(param_log_output_path_str, 'w') as usr_param_cfg:
                config.add_section('parameters')
                for param in list(hunter_param_dct.keys()):
                    config.set('parameters', str(param), str(hunter_param_dct[param]))
                config.write(usr_param_cfg)
                ready_to_run = True
        except IOError:
            error_log_lst.append('!! Failed to save parameter log files !!\n%s' % param_log_output_path_str)

        print(hunter_param_dct)

        error_log_lst = [_f for _f in error_log_lst if _f]

        if len(error_log_lst) > 0:
            print('Parameter error:', error_log_lst)
            error_log_lst.append('!!! Please check your settings !!!')
            self.ui.tab_a_statusrun_pte.appendPlainText('\n'.join(error_log_lst))
        else:
            if ready_to_run is True:
                self.ui.tab_a_runhunter_pb.setText(QtGui.QApplication.translate('MainWindow', 'Hunting ...', None,
                                                                                QtGui.QApplication.UnicodeUTF8))
                self.ui.tab_a_runhunter_pb.setEnabled(False)
                self.ui.tab_b_runbatch_pb.setEnabled(False)
                self.ui.tab_c_lmrun_pb.setEnabled(False)
                self.ui.tab_a_runhunter_pgb.setMinimum(0)
                self.ui.tab_a_runhunter_pgb.setMaximum(0)
                self.ui.tab_a_runhunter_pgb.show()
                self.single_worker.request_work(hunter_param_dct)

    def single_worker_info_update(self):

        back_info_str = self.single_worker.infoback()
        self.ui.tab_a_statusrun_pte.appendPlainText(back_info_str)
        # self.ui.tab_a_statusrun_pte.appendPlainText('\n')

    def batch_worker_on_finish(self):
        self.batch_thread.quit()
        print('!! batch_worker stopped !!')
        self.ui.tab_b_runbatch_pb.setText(QtGui.QApplication.translate('MainWindow',
                                                                       'Run batch mode identification >>>', None,
                                                                       QtGui.QApplication.UnicodeUTF8))
        self.ui.tab_b_runbatch_pgb.setMinimum(0)
        self.ui.tab_b_runbatch_pgb.setMaximum(100)
        self.ui.tab_b_runbatch_pgb.hide()
        self.ui.tab_b_runbatch_pb.setEnabled(True)
        self.ui.tab_a_runhunter_pb.setEnabled(True)
        self.ui.tab_c_lmrun_pb.setEnabled(True)

    def batch_worker_hunter(self):

        self.ui.tab_b_statusrun_pte.clear()
        ready_to_run = False

        loaded_cfg_files = str(self.ui.tab_b_infiles_pte.toPlainText())
        pre_loaded_cfg_lst = loaded_cfg_files.split('\n')

        # max_process = self.ui.tab_b_maxbatch_spb.value()
        sub_max_core = self.ui.tab_b_maxsubcore_spb.value()
        sub_max_ram = self.ui.tab_b_maxsubram_spb.value()

        loaded_cfg_lst = []
        for f in pre_loaded_cfg_lst:
            if len(f) > 4:
                loaded_cfg_lst.append(f)

        tot_num = len(loaded_cfg_lst)
        run_counter = 0

        os.chdir(self.lipidhunter_cwd)

        cfg_params_dct = {}
        self.ui.tab_b_statusrun_pte.appendPlainText('>>> Hunter started ... Total run number: %i \n\n' % tot_num)

        for _cfg in loaded_cfg_lst:

            hunter_param_dct, cfg_error = self.b_read_cfg(_cfg)
            if len(cfg_error) > 0:
                self.ui.tab_b_statusrun_pte.appendPlainText(str(cfg_error))
            if 'vendor' in list(hunter_param_dct.keys()):
                hunter_param_dct['batch_cfg_file'] = _cfg
                hunter_param_dct['core_number'] = sub_max_core
                hunter_param_dct['max_ram'] = sub_max_ram

                # capability issues to old batch mode files
                if 'lipid_class' in list(hunter_param_dct.keys()):
                    hunter_param_dct['lipid_type'] = hunter_param_dct['lipid_class']
                    ready_to_run = True
                elif 'lipid_type' in list(hunter_param_dct.keys()):
                    hunter_param_dct['lipid_class'] = hunter_param_dct['lipid_type']
                    self.ui.tab_b_statusrun_pte.appendPlainText(
                        '[WARNING] Please use "lipid_class" instead of "lipid_type" in your batch config file...\n')
                    ready_to_run = True
                else:
                    ready_to_run = False

                run_counter += 1
                cfg_params_dct[run_counter] = [hunter_param_dct, _cfg]

            else:
                run_counter += 1
                self.ui.tab_b_statusrun_pte.appendPlainText('!! Failed to read batch mode configure files: # %i / %i\n'
                                                            '%s\n!! Please check your settings '
                                                            '... skip this one ...\n' % (run_counter, tot_num, _cfg))
        if ready_to_run is True:

            self.ui.tab_b_runbatch_pb.setText(QtGui.QApplication.translate('MainWindow', 'Hunting in batch mode ...',
                                                                           None, QtGui.QApplication.UnicodeUTF8))
            self.ui.tab_b_runbatch_pb.setEnabled(False)
            self.ui.tab_a_runhunter_pb.setEnabled(False)
            self.ui.tab_c_lmrun_pb.setEnabled(False)
            self.ui.tab_b_runbatch_pgb.setMinimum(0)
            self.ui.tab_b_runbatch_pgb.setMaximum(0)
            self.ui.tab_b_runbatch_pgb.show()

            self.batch_worker.request_work(cfg_params_dct, tot_num)
        else:
            self.ui.tab_b_statusrun_pte.appendPlainText('!! Failed to read ALL batch mode configure files ...')

    def batch_worker_info_update(self):

        back_info_str = self.batch_worker.infoback()
        print('Got info: ', back_info_str)
        self.ui.tab_b_statusrun_pte.appendPlainText(back_info_str)

    def lm_worker_on_finish(self):
        self.lm_thread.quit()
        print('!! LipidMaster table export worker stopped !!')
        self.ui.tab_c_lmrun_pb.setText(QtGui.QApplication.translate('MainWindow',
                                                                    'Generate Lipid Master table >>>', None,
                                                                    QtGui.QApplication.UnicodeUTF8))
        self.ui.tab_b_runbatch_pb.setEnabled(True)
        self.ui.tab_a_runhunter_pb.setEnabled(True)
        self.ui.tab_c_lmrun_pb.setEnabled(True)
        self.ui.tab_c_runlm_pgb.setMinimum(0)
        self.ui.tab_c_runlm_pgb.setMaximum(100)
        self.ui.tab_c_runlm_pgb.hide()

    def lm_worker_hunter(self):

        self.ui.tab_c_lmstatus_pte.clear()

        _lipid_class_info = str(self.ui.tab_c_lipidclass_cmb.currentText())
        lipid_class_checker = re.compile(r'(.*)( [(])(\w{2,3})([)] )(.*)')
        lipid_class_match = lipid_class_checker.match(_lipid_class_info)

        lipid_class_chosen = False
        fa_list_chosen = False
        lm_export_chosen = False
        if lipid_class_match:
            lipid_class_info_lst = lipid_class_match.groups()
            usr_lipid_class = lipid_class_info_lst[2]
            usr_lipid_charge = lipid_class_info_lst[4]
            lipid_class_chosen = True
        else:
            usr_lipid_class = ''
            usr_lipid_charge = ''
            self.ui.tab_c_lmstatus_pte.appendPlainText('\n!! Please select a lipid class!!')

        usr_ms2_ppm = self.ui.tab_c_lmms2ppm_spb.value()

        fawhitelist_path_str = str(self.ui.tab_c_lmcalcfalist_le.text())
        lm_export_path_str = str(self.ui.tab_c_lmexport_le.text())
        if lm_export_path_str[-4:] != '.csv':
            lm_export_path_str += '.csv'
            self.ui.tab_c_lmexport_le.clear()
            self.ui.tab_c_lmexport_le.setText(lm_export_path_str)

        if os.path.isfile(fawhitelist_path_str):
            fa_list_chosen = True
        else:
            self.ui.tab_c_lmstatus_pte.appendPlainText('\n!! Please select a FA list!!')
        if len(lm_export_path_str) > 0:
            lm_export_chosen = True
        else:
            self.ui.tab_c_lmstatus_pte.appendPlainText('\n!! Please select a location to save the LipidMaster table!!')

        if lipid_class_chosen is True and fa_list_chosen is True and lm_export_chosen is True:

            os.chdir(self.lipidhunter_cwd)

            composer_param_dct = {'lipid_class': usr_lipid_class,
                                  'charge_mode': usr_lipid_charge,
                                  'exact_position': 'FALSE',
                                  'ms2ppm': usr_ms2_ppm,
                                  'fa_whitelist': fawhitelist_path_str,
                                  'export_path': lm_export_path_str,
                                  }
            print(composer_param_dct)
            self.ui.tab_c_lmstatus_pte.appendPlainText('>>> Start to generate LipidMaster table...')

            self.ui.tab_c_lmrun_pb.setText(QtGui.QApplication.translate('MainWindow', '... Generating ...',
                                                                        None, QtGui.QApplication.UnicodeUTF8))
            self.ui.tab_b_runbatch_pb.setEnabled(False)
            self.ui.tab_a_runhunter_pb.setEnabled(False)
            self.ui.tab_c_lmrun_pb.setEnabled(False)

            self.ui.tab_c_runlm_pgb.setMinimum(0)
            self.ui.tab_c_runlm_pgb.setMaximum(0)
            self.ui.tab_c_runlm_pgb.show()

            self.lm_worker.request_work(composer_param_dct)
        else:
            self.ui.tab_c_lmstatus_pte.appendPlainText('!! Please check your settings and try again !!\n')

    def lm_worker_info_update(self):

        back_info_str = self.lm_worker.infoback()
        print('Got info: ', back_info_str)
        self.ui.tab_c_lmstatus_pte.appendPlainText(back_info_str)


class SingleWorker(QtCore.QObject):
    workRequested = QtCore.Signal()
    finished = QtCore.Signal()
    info_update = QtCore.Signal(str)

    def __init__(self, parent=None):
        super(SingleWorker, self).__init__(parent)
        self.params_dct = {}
        self.run_count = 0
        self.total_count = 0
        self.info_str = ''

    def request_work(self, params_dct):
        """
        Get parameters from Main window
        :param(dict) params_dct: a dict contains all parameters from UI
        """
        self.workRequested.emit()
        self.params_dct = params_dct

    def infoback(self):

        return self.info_str

    def run_hunter(self):

        log_lst = []
        print('>>> Hunter single worker started ...')
        time.sleep(1)  # Wait for 3 sec to avoid overwriting the self.info_str
        self.info_str = '>>> Hunter started ... Please wait ...\n'
        self.infoback()
        self.info_update.emit(self.info_str)

        hunter_time, log_lst, output_df2 = huntlipids(self.params_dct, error_lst=log_lst)

        # try:
        #     hunter_time, log_lst, output_df2 = huntlipids(self.params_dct, error_lst=log_lst)
        # except Exception as _e:
        #     print(_e)
        #     hunter_time = False
        #     log_lst = False
        #     export_df = False
        #     time.sleep(1)
        #     self.info_str = '!! Sorry, an error has occurred, please check your settings !!'
        #     self.infoback()
        #     self.info_update.emit(self.info_str)
        #     time.sleep(1)
        #     self.finished.emit()

        err_info = ''
        if isinstance(log_lst, list):
            err_info = '\n'.join(log_lst)
        else:
            pass

        if hunter_time is not False:

            if isinstance(hunter_time, float):
                time.sleep(1)
                self.info_str = '%s\n>>> >>> >>> FINISHED in %.3f Sec <<< <<< <<<\n' % (err_info, hunter_time)
                self.infoback()
                self.info_update.emit(self.info_str)
            else:
                if isinstance(hunter_time, str):
                    time.sleep(1)
                    self.info_str = '%s\n>>> >>> >>> FINISHED in %s Sec <<< <<< <<<\n' % (err_info, hunter_time)
                    self.infoback()
                    self.info_update.emit(self.info_str)
                else:
                    time.sleep(1)
                    self.info_str = '%s\n!! Sorry, an error has occurred, please check your settings !!\n\n' % err_info
                    self.infoback()
                    self.info_update.emit(self.info_str)

        else:
            time.sleep(1)
            self.info_str = '%s\n!! Sorry, an error has occurred, please check your settings !!\n\n' % err_info
            self.infoback()
            self.info_update.emit(self.info_str)
            time.sleep(1)
            self.finished.emit()

        time.sleep(1)
        self.finished.emit()


class BatchWorker(QtCore.QObject):
    workRequested = QtCore.Signal()
    finished = QtCore.Signal()
    info_update = QtCore.Signal(str)

    def __init__(self, parent=None):
        super(BatchWorker, self).__init__(parent)
        self.cfg_params_dct = {}
        self.run_count = 0
        self.total_count = 0
        self.info_str = ''

    def request_work(self, cfg_params_dct, total_count):
        """
        Get parameters from Main window
        :param(dict) cfg_params_dct: a dict contains all configurations from all batch files
        :param(int) total_count: total number of config.txt need to run
        """
        self.workRequested.emit()
        self.cfg_params_dct = cfg_params_dct
        self.total_count = total_count

    def infoback(self):

        return self.info_str

    def run_hunter(self):

        t_start = time.time()

        cfg_key_lst = list(self.cfg_params_dct.keys())

        for cfg_idx in cfg_key_lst:
            log_lst = []
            ready_to_run = False
            time.sleep(3)  # Wait for 3 sec to avoid overwriting the self.info_str
            _param_dct = self.cfg_params_dct[cfg_idx][0]
            _cfg_path = self.cfg_params_dct[cfg_idx][1]
            self.run_count = int(cfg_idx)
            print('>>> Hunter batch_worker started ...')
            time.sleep(1)
            self.info_str = 'Start processing file: # %i / %i\n%s\n' % (cfg_idx, self.total_count, _cfg_path)
            self.infoback()
            self.info_update.emit(self.info_str)

            start_time_str = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
            _param_dct['hunter_start_time'] = start_time_str

            if not os.path.isdir(_param_dct['img_output_folder_str']):
                os.makedirs(_param_dct['img_output_folder_str'])
                self.info_str = 'Output folder created...\n\n'
                self.infoback()
                self.info_update.emit(self.info_str)
                print('Output folder created...')
            output_folder_path = os.path.abspath(_param_dct['img_output_folder_str'])
            log_file_name = 'LipidHunter_Params-Log_%s.txt' % _param_dct['hunter_start_time']
            param_log_output_path_str = os.path.join(output_folder_path, log_file_name)

            try:
                config = configparser.ConfigParser()
                with open(param_log_output_path_str, 'w') as usr_param_cfg:
                    config.add_section('parameters')
                    for param in list(_param_dct.keys()):
                        config.set('parameters', str(param), str(_param_dct[param]))
                    config.write(usr_param_cfg)
                try:
                    os.chdir(_param_dct['hunter_folder'])
                    hunter_py_path = os.path.join(_param_dct['hunter_folder'], 'LipidHunter.py')
                    hunter_exe_path = os.path.join(_param_dct['hunter_folder'], 'LipidHunter.exe')
                    if os.path.isfile(hunter_py_path):
                        print('>>> Running LipidHunter source code version ...')
                        ready_to_run = True
                    elif os.path.isfile(hunter_exe_path):
                        print('>>> Running LipidHunter .exe version ...')
                        ready_to_run = True
                    else:
                        ready_to_run = False
                        time.sleep(1)
                        self.info_str = '!! LipidHunter folder path in configuration is not correct!!\n'
                        self.infoback()
                        self.info_update.emit(self.info_str)
                except IOError:
                    ready_to_run = False
                    time.sleep(1)
                    self.info_str = '!! LipidHunter folder path in configuration is not correct!!\n'
                    self.infoback()
                    self.info_update.emit(self.info_str)
            except IOError:
                ready_to_run = False
                time.sleep(1)
                self.info_str = '!! Failed to save parameter log files ...\n'
                self.infoback()
                self.info_update.emit(self.info_str)

            if ready_to_run is True:
                hunter_time, log_lst, output_df2 = huntlipids(_param_dct, error_lst=log_lst)
                # try:
                #     hunter_time, log_lst, output_df2 = huntlipids(_param_dct, error_lst=log_lst)
                # except Exception as _e:
                #     print(_e)
                #     hunter_time = False
                #     log_lst = False
                #     export_df = False

                err_info = ''
                if isinstance(log_lst, list):
                    err_info = '\n'.join(log_lst)
                else:
                    pass

                if hunter_time is not False:

                    print('Hunter finished successfully ...')
                    if isinstance(hunter_time, float):
                        run_time = '%.3f' % hunter_time
                        time.sleep(1)
                        self.info_str = ('%s\n>>> FINISHED with file # %i / %i in %s Sec ...\n\n'
                                         % (err_info, cfg_idx, self.total_count, run_time))
                        self.infoback()
                        self.info_update.emit(self.info_str)

                    elif isinstance(hunter_time, str):
                        time.sleep(1)
                        self.info_str = ('%s\n>>>FINISHED with file # %i / %i in %s Sec ...\n\n'
                                         % (err_info, cfg_idx, self.total_count, hunter_time))
                        self.infoback()
                        self.info_update.emit(self.info_str)

                    else:
                        err_info = ''
                        for err in log_lst:
                            err_info += err
                        time.sleep(1)
                        self.info_str = ('%s\n!! Failed to process batch mode configuration file ... '
                                         'skip this one ...\n\n' % err_info)
                        self.infoback()
                        self.info_update.emit(self.info_str)

                else:
                    time.sleep(1)
                    self.info_str = ('%s\n!! Failed to process batch mode configuration file ... skip this one ...\n\n'
                                     % err_info)
                    self.infoback()
                    self.info_update.emit(self.info_str)

            else:
                time.sleep(1)
                self.info_str = ('!! Failed to save parameter log files ...\n'
                                 '!! Failed to process batch mode configuration file ... skip this one ...\n\n')
                self.infoback()
                self.info_update.emit(self.info_str)

        t_end = time.time() - t_start
        time.sleep(3)
        self.info_str = '\n>>> ALL FINISHED in %.3f Sec <<<' % t_end
        self.infoback()
        self.info_update.emit(self.info_str)
        self.finished.emit()


class LMWorker(QtCore.QObject):
    workRequested = QtCore.Signal()
    finished = QtCore.Signal()
    info_update = QtCore.Signal(str)

    def __init__(self, parent=None):
        super(LMWorker, self).__init__(parent)
        self.lm_params_dct = {}
        self.info_str = ''

    def request_work(self, lm_params_dct):
        """
        Get parameters from Main window
        :param(dict) lm_params_dct: a dict contains all params for LipidMaster table generator
        """
        self.workRequested.emit()
        self.lm_params_dct = lm_params_dct

    def infoback(self):
        return self.info_str

    def run_lm_generator(self):

        t_start = time.time()

        lipidcomposer = LipidComposer()
        usr_lipid_master_df = lipidcomposer.compose_lipid(param_dct=self.lm_params_dct,
                                                          ms2_ppm=self.lm_params_dct['ms2ppm'])

        if isinstance(usr_lipid_master_df, pd.DataFrame):
            if not usr_lipid_master_df.empty:
                self.info_str = ('==> Number of predicted lipids (discrete form): %i' % usr_lipid_master_df.shape[0])
                self.infoback()
                self.info_update.emit(self.info_str)

                abs_export_path = os.path.abspath(self.lm_params_dct['export_path'])
                abs_export_folder = os.path.dirname(abs_export_path)
                if os.path.isdir(abs_export_folder):
                    pass
                else:
                    os.mkdir(os.path.abspath(abs_export_folder))
                try:
                    usr_lipid_master_df.to_csv(abs_export_path)
                    self.info_str = ('==> --> Lipid Master table Saved as: \n %s' % abs_export_path)
                    self.infoback()
                    self.info_update.emit(self.info_str)

                except Exception as _err:
                    print(_err)
                    self.info_str = '\n %s \n' % _err
                    self.infoback()
                    self.info_update.emit(self.info_str)

            else:
                print('[ERROR] Failed to generate LipidMaster Table...\n'
                      '... ... Check if Lipid Class and FA are marked in FA whitelist...')
                self.info_str = ('[ERROR] Failed to generate LipidMaster Table...\n'
                                 ' ... ... Please check if Lipid Class and FA are marked in FA whitelist...')
                self.infoback()
                self.info_update.emit(self.info_str)
        else:
            print('[ERROR] Failed to generate LipidMaster Table...\n'
                  '... ... Check if Lipid Class and FA are marked in FA whitelist...')
            self.info_str = ('[ERROR] Failed to generate LipidMaster Table...\n'
                             '... ... Please check if Lipid Class and FA are marked in FA whitelist...')
            self.infoback()
            self.info_update.emit(self.info_str)

        t_end = time.time() - t_start

        time.sleep(0.25)
        self.lm_params_dct = {}
        self.info_str = '\n>>> FINISHED in %.3f Sec <<<' % t_end
        self.infoback()
        self.info_update.emit(self.info_str)
        self.finished.emit()


if __name__ == '__main__':
    import sys

    multiprocessing.freeze_support()
    app = QtGui.QApplication(sys.argv)
    window = LipidHunterMain()
    window.show()
    sys.exit(app.exec_())
