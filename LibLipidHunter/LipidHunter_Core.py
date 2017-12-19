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
# try:  # python3
#     import configparser
# except NameError:  # python2
#     import ConfigParser as configparser

from __future__ import print_function
import ConfigParser as configparser
import glob
import os
import re
import time

import pandas as pd
from PySide import QtCore, QtGui

from LibLipidHunter import ExtractorMZML
from LibLipidHunter.LinkerMZML import hunt_link
from LibLipidHunter.LipidHunter_UI import Ui_MainWindow
from LibLipidHunter.SpectraHunter import huntlipids


class LipidHunterCore(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None, cwd=None):
        QtGui.QMainWindow.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        # set version
        version_date = r'19, December, 2017'
        version_html = (r'<html><head/><body><p><span style=" font-weight:600;">'
                        r'LipidHunter Beta released date: {version_date}'
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

        # slots for tab a
        QtCore.QObject.connect(self.ui.tab_a_loadxlsxpath_pb, QtCore.SIGNAL("clicked()"), self.a_load_xlsx)
        QtCore.QObject.connect(self.ui.tab_a_launchgen_pb, QtCore.SIGNAL("clicked()"), self.a_go_generator)
        QtCore.QObject.connect(self.ui.tab_a_mzml_pb, QtCore.SIGNAL("clicked()"), self.a_load_mzml)
        QtCore.QObject.connect(self.ui.tab_a_saveimgfolder_pb, QtCore.SIGNAL("clicked()"), self.a_save_img2folder)
        QtCore.QObject.connect(self.ui.tab_a_sumxlsxpath_pb, QtCore.SIGNAL("clicked()"), self.a_save_output)
        QtCore.QObject.connect(self.ui.tab_a_runhunter_pb, QtCore.SIGNAL("clicked()"), self.a_run_hunter)
        # # slots for tab c
        QtCore.QObject.connect(self.ui.tab_c_fawhitelist_pb, QtCore.SIGNAL("clicked()"), self.c_load_fawhitelist)
        QtCore.QObject.connect(self.ui.tab_c_hgcfg_pb, QtCore.SIGNAL("clicked()"), self.c_load_hgcfg)
        QtCore.QObject.connect(self.ui.tab_c_scorecfg_pb, QtCore.SIGNAL("clicked()"), self.c_load_scorecfg)
        QtCore.QObject.connect(self.ui.tab_c_savesettings_pb, QtCore.SIGNAL("clicked()"), self.c_set_default_cfg)

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
            if 'fa_white_list_cfg' in options:
                self.ui.tab_c_fawhitelist_le.setText(config.get(user_cfg, 'fa_white_list_cfg'))
            if 'lipid_specific_cfg' in options:
                self.ui.tab_c_hgcfg_le.setText(config.get(user_cfg, 'lipid_specific_cfg'))
            if 'score_cfg' in options:
                self.ui.tab_c_scorecfg_le.setText(config.get(user_cfg, 'score_cfg'))
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
        if folder is not u'':
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
    def check_path(usr_path, info_str):
        try:
            if os.path.isfile(usr_path):
                error_log = ''
            else:
                error_log = '!! Failed to load {_file} !!'.format(_file=info_str)
        except IOError:
            error_log = '!! Failed to load {_file} !!'.format(_file=info_str)

        return error_log

    def a_load_xlsx(self):
        file_info_str = u'FA white list files (*.xlsx *.XLSX)'
        self.open_file(file_info_str, self.ui.tab_a_loadxlsxpath_le)

    def a_go_generator(self):
        self.ui.tabframe.setCurrentIndex(3)

    def a_load_mzml(self):
        file_info_str = u'mzML spectra files (*.mzML *.mzml)'
        self.open_file(file_info_str, self.ui.tab_a_mzml_le)

    def a_save_img2folder(self):
        a_save_img2folder_str = QtGui.QFileDialog.getExistingDirectory()
        self.ui.tab_a_saveimgfolder_le.clear()
        self.ui.tab_a_saveimgfolder_le.setText(a_save_img2folder_str)

    def a_save_output(self):
        a_save_output_path = QtGui.QFileDialog.getSaveFileName(caption=u'Save file', filter=u'.xlsx')
        self.ui.tab_a_savexlsxpath_le.clear()
        a_save_output_str = os.path.abspath(a_save_output_path[0])
        self.ui.tab_a_savexlsxpath_le.setText(a_save_output_str)

    def a_run_hunter(self):

        error_log_lst = []

        usr_vendor_str = self.ui.vendor_cmb.currentText().lower()
        vendors_dct = {'therm': 'thermo', 'water': 'waters', 'sciex': 'sciex', 'agile': 'agilent'}
        if usr_vendor_str[0:5] in vendors_dct.keys():
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

        _pl_class_info = str(self.ui.tab_a_lipidclass_cmb.currentText())

        pl_class_checker = re.compile(r'(.*)( [\(])(\w{2,3})([\)] )(.*)')

        pl_class_match = pl_class_checker.match(_pl_class_info)

        if pl_class_match:
            pl_class_info_lst = pl_class_match.groups()
            _pl_class = pl_class_info_lst[2]
            _pl_charge = pl_class_info_lst[4]
        else:
            _pl_class = ''
            _pl_charge = ''
            error_log_lst.append('!! Please select a lipid class!!')

        fawhitelist_path_str = str(self.ui.tab_a_loadxlsxpath_le.text())
        mzml_path_str = str(self.ui.tab_a_mzml_le.text())
        img_output_folder_str = str(self.ui.tab_a_saveimgfolder_le.text())
        xlsx_output_path_str = str(self.ui.tab_a_savexlsxpath_le.text())

        rt_start = self.ui.tab_a_rtstart_dspb.value()
        rt_end = self.ui.tab_a_rtend_dspb.value()
        mz_start = self.ui.tab_a_mzstart_dspb.value()
        mz_end = self.ui.tab_a_mzend_dspb.value()
        dda_top = self.ui.tab_a_dda_spb.value()
        ms_th = self.ui.tab_a_msthreshold_spb.value()
        ms2_th = self.ui.tab_a_ms2threshold_spb.value()
        ms_ppm = self.ui.tab_a_msppm_spb.value()
        ms2_ppm = self.ui.tab_a_ms2ppm_spb.value()
        hg_th = self.ui.tab_a_hgthreshold_spb.value()
        hg_ppm = self.ui.tab_a_hgppm_spb.value()
        score_filter = self.ui.tab_a_score_spb.value()
        isotope_score_filter = self.ui.tab_a_isotopescore_spb.value()
        ms2_info_threshold = self.ui.tab_a_ms2infoth_dspb.value() * 0.01
        hgms2_info_threshold = self.ui.tab_a_ms2hginfoth_dspb.value() * 0.01

        fa_white_list_cfg = self.ui.tab_c_fawhitelist_le.text()
        lipid_specific_cfg = self.ui.tab_c_hgcfg_le.text()
        score_cfg = self.ui.tab_c_scorecfg_le.text()

        print('Vendor mode = %s, Experiment mode = %s' % (usr_vendor, usr_exp_mode))
        print('Hunter started!')

        error_log_lst.append(self.check_path(fawhitelist_path_str, 'FA whitelist'))
        error_log_lst.append(self.check_path(mzml_path_str, 'mzML spectra'))
        error_log_lst.append(self.check_path(fa_white_list_cfg, 'FA whitelist'))
        error_log_lst.append(self.check_path(lipid_specific_cfg, 'configuration for Phospholipids'))
        error_log_lst.append(self.check_path(score_cfg, 'W_frag score configuration'))

        usr_score_mode = self.ui.tab_c_scoremode_cmb.currentIndex()
        if usr_score_mode == 0:
            print(self.ui.tab_c_scoremode_cmb.currentText())
            rank_score = True
        else:
            print(self.ui.tab_c_scoremode_cmb.currentText())
            rank_score = False

        usr_isotope_score_mode = self.ui.tab_c_isotopescoremode_cmb.currentIndex()
        if usr_isotope_score_mode == 0:
            print(self.ui.tab_c_isotopescoremode_cmb.currentText())
            fast_isotope = False
        else:
            print(self.ui.tab_c_isotopescoremode_cmb.currentText())
            fast_isotope = True

        start_time_str = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())

        hunter_param_dct = {'fawhitelist_path_str': fawhitelist_path_str, 'mzml_path_str': mzml_path_str,
                            'img_output_folder_str': img_output_folder_str,
                            'xlsx_output_path_str': xlsx_output_path_str, 'rt_start': rt_start, 'rt_end': rt_end,
                            'mz_start': mz_start, 'mz_end': mz_end, 'dda_top': dda_top, 'ms_th': ms_th,
                            'ms2_th': ms2_th, 'ms_ppm': ms_ppm, 'ms2_ppm': ms2_ppm, 'hg_th': hg_th, 'hg_ppm': hg_ppm,
                            'score_filter': score_filter, 'isotope_score_filter': isotope_score_filter,
                            'lipid_type': _pl_class, 'charge_mode': _pl_charge, 'fa_white_list_cfg': fa_white_list_cfg,
                            'lipid_specific_cfg': lipid_specific_cfg, 'score_cfg': score_cfg, 'vendor': usr_vendor,
                            'ms2_infopeak_threshold': ms2_info_threshold,
                            'ms2_hginfopeak_threshold': hgms2_info_threshold,
                            'rank_score': rank_score, 'fast_isotope': fast_isotope,
                            'hunter_folder': self.lipidhunter_cwd,
                            'hunter_start_time': start_time_str, 'experiment_mode': usr_exp_mode}

        param_log_output_path_str = (str(self.ui.tab_a_saveimgfolder_le.text()) +
                                     '/LipidHunter_Params-Log_%s.txt' % start_time_str
                                     )

        try:
            config = configparser.ConfigParser()
            with open(param_log_output_path_str, 'w') as usr_param_cfg:
                config.add_section('parameters')
                for param in hunter_param_dct.keys():
                    config.set('parameters', param, hunter_param_dct[param])
                config.write(usr_param_cfg)

        except IOError:
            error_log_lst.append('!! Failed to save parameter log files !!')

        print(hunter_param_dct)

        error_log_lst = filter(None, error_log_lst)
        print(error_log_lst)
        if len(error_log_lst) > 0:
            error_log_lst.append('!!! Please check your settings !!!')
            self.ui.tab_a_statusrun_pte.appendPlainText('\n'.join(error_log_lst) + '\n')
        else:
            # default output code
            try:
                tot_run_time = huntlipids(hunter_param_dct)

            except:
                tot_run_time = '!! Sorry, an error has occurred, please check your settings !!'

            if isinstance(tot_run_time, float):
                self.ui.tab_a_statusrun_pte.insertPlainText('%.2f Sec\n' % tot_run_time)
                self.ui.tab_a_statusrun_pte.insertPlainText('>>> >>> >>> FINISHED <<< <<< <<<')

            else:
                if isinstance(tot_run_time, str):
                    self.ui.tab_a_statusrun_pte.appendPlainText(tot_run_time)

                else:
                    self.ui.tab_a_statusrun_pte.appendPlainText('!! Sorry, an error has occurred, '
                                                                'please check your settings !!')

            # # for debug only
            # tot_run_time = huntlipids(hunter_param_dct)
            # self.ui.tab_e_statusrun_pte.insertPlainText('%.2f Sec\n' % tot_run_time)
            # self.ui.tab_e_statusrun_pte.insertPlainText('>>> >>> >>> FINISHED <<< <<< <<<')

    def c_set_default_cfg(self):
        config = configparser.ConfigParser()
        with open('config.ini', 'w') as default_cfg:
            config.add_section('settings')
            config.set('settings', 'fa_white_list_cfg', self.ui.tab_c_fawhitelist_le.text())
            config.set('settings', 'lipid_specific_cfg', self.ui.tab_c_hgcfg_le.text())
            config.set('settings', 'score_cfg', self.ui.tab_c_scorecfg_le.text())

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

    def c_load_fawhitelist(self):
        file_info_str = u'FA white list files (*.xlsx *.XLSX)'
        self.open_file(file_info_str, self.ui.tab_c_fawhitelist_le)

    def c_load_hgcfg(self):
        file_info_str = u'MS Excel files (*.xlsx *.XLSX)'
        self.open_file(file_info_str, self.ui.tab_c_hgcfg_le)

    def c_load_scorecfg(self):
        file_info_str = u'MS Excel files (*.xlsx *.XLSX)'
        self.open_file(file_info_str, self.ui.tab_c_scorecfg_le)


if __name__ == '__main__':
    import sys

    app = QtGui.QApplication(sys.argv)
    window = LipidHunterCore()
    window.show()
    sys.exit(app.exec_())
