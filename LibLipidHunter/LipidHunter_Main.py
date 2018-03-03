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

from __future__ import division
from __future__ import print_function

import glob
import multiprocessing
import multiprocessing.pool
import os
import re
import time

from PySide import QtCore, QtGui
from six.moves import configparser

try:
    # import configparser
    from LibLipidHunter.LipidHunter_UI import Ui_MainWindow
    from LibLipidHunter.Hunter_Core import huntlipids
except ImportError:  # for python 2.7.14
    # import ConfigParser as configparser
    from LipidHunter_UI import Ui_MainWindow
    from Hunter_Core import huntlipids


class LipidHunterMain(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None, cwd=None):
        QtGui.QMainWindow.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        # set version
        version_date = r'02, March, 2018'
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
        self.a_max_ms()

        # disable multi_mode in batch run
        self.ui.tab_b_mutlimode_cmb.hide()
        self.ui.tab_b_maxbatch_lb.hide()
        self.ui.tab_b_maxbatch_spb.setValue(1)
        self.ui.tab_b_maxbatch_spb.hide()
        self.ui.tab_b_maxsubcore_spb.setValue(3)
        self.ui.tab_b_maxsubram_spb.setValue(5)

        # slots for tab a
        QtCore.QObject.connect(self.ui.tab_a_loadxlsxpath_pb, QtCore.SIGNAL("clicked()"), self.a_load_xlsx)
        QtCore.QObject.connect(self.ui.tab_a_launchgen_pb, QtCore.SIGNAL("clicked()"), self.a_go_generator)
        QtCore.QObject.connect(self.ui.tab_a_mzml_pb, QtCore.SIGNAL("clicked()"), self.a_load_mzml)
        QtCore.QObject.connect(self.ui.tab_a_saveimgfolder_pb, QtCore.SIGNAL("clicked()"), self.a_save_img2folder)
        QtCore.QObject.connect(self.ui.tab_a_msmax_chb, QtCore.SIGNAL("clicked()"), self.a_max_ms)
        QtCore.QObject.connect(self.ui.tab_a_sumxlsxpath_pb, QtCore.SIGNAL("clicked()"), self.a_save_output)
        QtCore.QObject.connect(self.ui.tab_a_runhunter_pb, QtCore.SIGNAL("clicked()"), self.a_run_hunter)
        QtCore.QObject.connect(self.ui.tab_a_cfgpath_pb, QtCore.SIGNAL("clicked()"), self.a_save_cfg)
        QtCore.QObject.connect(self.ui.tab_a_gencfg_pb, QtCore.SIGNAL("clicked()"), self.a_create_cfg)
        # # slots for tab b
        self.ui.tab_b_mutlimode_cmb.currentIndexChanged['QString'].connect(self.b_set_multi_mode)
        QtCore.QObject.connect(self.ui.tab_b_addcfg_pb, QtCore.SIGNAL("clicked()"), self.b_load_batchcfg)
        QtCore.QObject.connect(self.ui.tab_b_addcfgfolder_pb, QtCore.SIGNAL("clicked()"), self.b_load_batchcfgfolder)
        QtCore.QObject.connect(self.ui.tab_b_clearall_pb, QtCore.SIGNAL("clicked()"), self.ui.tab_b_infiles_pte.clear)
        QtCore.QObject.connect(self.ui.tab_b_runbatch_pb, QtCore.SIGNAL("clicked()"), self.b_run_batchmode)
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
                self.ui.tab_a_loadxlsxpath_le.setText(config.get(user_cfg, 'fa_white_list_cfg'))
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
        try:
            if os.path.isfile(usr_path):
                error_log = ''
            else:
                error_log = '!! Failed to load {_file} !!'.format(_file=info_str)
        except IOError:
            error_log = '!! Failed to load {_file} !!'.format(_file=info_str)

        return error_log

    @staticmethod
    def check_folder(usr_path, info_str):
        try:
            if os.path.isdir(usr_path):
                error_log = ''
            else:
                os.makedirs(usr_path)
                print('Folder created... %s' % usr_path)
                error_log = ''
        except IOError:
            error_log = '!! Failed to open folder {_file} !!'.format(_file=info_str)

        return error_log

    def a_load_xlsx(self):
        file_info_str = 'FA white list files (*.xlsx *.XLSX)'
        self.open_file(file_info_str, self.ui.tab_a_loadxlsxpath_le)

    def a_go_generator(self):
        self.ui.tabframe.setCurrentIndex(3)

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
        a_save_output_str = os.path.abspath(a_save_output_path[0])
        self.ui.tab_a_savexlsxpath_le.setText(a_save_output_str)

    def a_max_ms(self):
        if self.ui.tab_a_msmax_chb.isChecked():
            self.ui.tab_a_msmax_spb.show()
            self.ui.tab_a_msmax_spb.setValue(100 * self.ui.tab_a_msthreshold_spb.value())
        else:
            self.ui.tab_a_msmax_spb.hide()
            self.ui.tab_a_msmax_spb.setValue(0)

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

        _pl_class_info = str(self.ui.tab_a_lipidclass_cmb.currentText())

        pl_class_checker = re.compile(r'(.*)( [(])(\w{2,3})([)] )(.*)')

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
        img_output_folder_str = str(self.ui.tab_a_saveimgfolder_le.text()).strip(r'\/')
        xlsx_output_path_str = str(self.ui.tab_a_savexlsxpath_le.text())

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
        hg_th = self.ui.tab_a_hgthreshold_spb.value()
        hg_ppm = self.ui.tab_a_hgppm_spb.value()
        ms2_info_threshold = self.ui.tab_a_ms2infoth_dspb.value() * 0.01
        hgms2_info_threshold = self.ui.tab_a_ms2hginfoth_dspb.value() * 0.01

        ms_max = 0
        if self.ui.tab_a_msmax_chb.isChecked() and self.ui.tab_a_msmax_spb.value() > ms_th + 1:
            ms_max = self.ui.tab_a_msmax_spb.value()
        lipid_specific_cfg = self.ui.tab_c_hgcfg_le.text()
        score_cfg = self.ui.tab_c_scorecfg_le.text()

        core_num = self.ui.tab_c_cores_spb.value()
        max_ram = self.ui.tab_c_ram_spb.value()
        img_typ = self.ui.tab_c_imagetype_cmb.currentText()[1:]
        img_dpi = self.ui.tab_c_dpi_spb.value()

        error_log_lst.append(self.check_file(fawhitelist_path_str, 'FA whitelist'))
        error_log_lst.append(self.check_file(mzml_path_str, 'mzML spectra'))
        error_log_lst.append(self.check_file(lipid_specific_cfg, 'configuration for Phospholipids'))
        error_log_lst.append(self.check_file(score_cfg, 'W_frag score configuration'))
        error_log_lst.append(self.check_folder(img_output_folder_str, 'Output folder'))

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
            tag_all_sn = False
        else:
            tag_all_sn = True

        start_time_str = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())

        score_filter = rank_score_filter

        hunter_param_dct = {'fawhitelist_path_str': fawhitelist_path_str, 'mzml_path_str': mzml_path_str,
                            'img_output_folder_str': img_output_folder_str,
                            'xlsx_output_path_str': xlsx_output_path_str, 'rt_start': rt_start, 'rt_end': rt_end,
                            'mz_start': mz_start, 'mz_end': mz_end, 'dda_top': dda_top, 'ms_th': ms_th,
                            'ms2_th': ms2_th, 'ms_ppm': ms_ppm, 'ms2_ppm': ms2_ppm, 'hg_th': hg_th, 'hg_ppm': hg_ppm,
                            'rank_score_filter': rank_score_filter, 'isotope_score_filter': isotope_score_filter,
                            'score_filter': score_filter,
                            'lipid_type': _pl_class, 'charge_mode': _pl_charge,
                            'lipid_specific_cfg': lipid_specific_cfg, 'score_cfg': score_cfg, 'vendor': usr_vendor,
                            'ms2_infopeak_threshold': ms2_info_threshold,
                            'ms2_hginfopeak_threshold': hgms2_info_threshold,
                            'rank_score': rank_score, 'fast_isotope': fast_isotope,
                            'hunter_folder': self.lipidhunter_cwd,
                            'hunter_start_time': start_time_str, 'experiment_mode': usr_exp_mode,
                            'core_number': core_num, 'max_ram': max_ram, 'tag_all_sn': tag_all_sn,
                            'img_type': img_typ, 'img_dpi': img_dpi, 'ms_max': ms_max, 'pr_window': pr_window}

        return hunter_param_dct, error_log_lst

    def a_run_hunter(self):

        self.ui.tab_a_statusrun_pte.clear()
        self.ui.tab_a_statusrun_pte.setPlainText('')

        hunter_param_dct, error_log_lst = self.a_get_params()

        print('Vendor mode = %s, Experiment mode = %s' % (hunter_param_dct['vendor'],
                                                          hunter_param_dct['experiment_mode']))
        print('Isotope score mode = %s' % (hunter_param_dct['fast_isotope']))
        print('Rankscore mode = %s' % (hunter_param_dct['rank_score']))
        print('Hunter started!')

        output_folder_path = str(self.ui.tab_a_saveimgfolder_le.text()).strip(r'\/')

        if os.path.isdir(output_folder_path):
            print('Output folder path... %s' % output_folder_path)
        else:
            try:
                os.mkdir(output_folder_path)
                print('Output folder created... %s' % output_folder_path)
            except IOError:
                error_log_lst.append('!! Failed to create output folder !!')

        param_log_output_path_str = (output_folder_path + '/LipidHunter_Params-Log_%s.txt'
                                     % hunter_param_dct['hunter_start_time'])

        try:
            config = configparser.ConfigParser()
            with open(param_log_output_path_str, 'w') as usr_param_cfg:
                config.add_section('parameters')
                for param in list(hunter_param_dct.keys()):
                    config.set('parameters', str(param), str(hunter_param_dct[param]))
                config.write(usr_param_cfg)

        except IOError:
            error_log_lst.append('!! Failed to save parameter log files !!')

        print(hunter_param_dct)

        error_log_lst = [_f for _f in error_log_lst if _f]

        if len(error_log_lst) > 0:
            print('Parameter error:', error_log_lst)
            error_log_lst.append('!!! Please check your settings !!!')
            self.ui.tab_a_statusrun_pte.appendPlainText('\n'.join(error_log_lst) + '\n')
        else:

            # # for debug only
            # tot_run_time = huntlipids(hunter_param_dct)
            # self.ui.tab_a_statusrun_pte.insertPlainText('%.2f Sec\n' % tot_run_time)
            # self.ui.tab_a_statusrun_pte.insertPlainText('>>> >>> >>> FINISHED <<< <<< <<<')

            # default output code
            try:
                tot_run_time, error_log_lst, export_df = huntlipids(hunter_param_dct, error_log_lst)

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
                    if len(error_log_lst) > 0:
                        for err in error_log_lst:
                            self.ui.tab_a_statusrun_pte.appendPlainText(str(err) + '\n')

    def a_save_cfg(self):
        a_save_cfg_path = QtGui.QFileDialog.getSaveFileName(caption='Save file', filter='.txt')
        self.ui.tab_a_cfgpath_le.clear()
        a_save_cfg_str = os.path.abspath(a_save_cfg_path[0])
        self.ui.tab_a_cfgpath_le.setText(a_save_cfg_str)

    def a_create_cfg(self):
        hunter_param_dct, error_log_lst = self.a_get_params()
        error_log_lst = [_f for _f in error_log_lst if _f]

        if len(error_log_lst) > 0:
            print('Parameter error:', error_log_lst)
            error_log_lst.append('!!! Please check your settings !!!')
            self.ui.tab_a_gencfg_pte.appendPlainText('\n'.join(error_log_lst) + '\n')
        else:
            param_cfg_path_str = str(self.ui.tab_a_cfgpath_le.text())
            param_cfg_directory = os.path.dirname(param_cfg_path_str)
            if not os.path.exists(param_cfg_directory):
                os.makedirs(param_cfg_directory)
            config = configparser.ConfigParser()
            with open(param_cfg_path_str, 'w') as usr_param_cfg:
                config.add_section('parameters')
                for param in list(hunter_param_dct.keys()):
                    config.set('parameters', str(param), str(hunter_param_dct[param]))
                config.write(usr_param_cfg)
                self.ui.tab_a_gencfg_pte.insertPlainText('>>> Configuration saved as:')
                self.ui.tab_a_gencfg_pte.insertPlainText(param_cfg_path_str)
                self.ui.tab_a_gencfg_pte.insertPlainText('\n')

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

        i_type_key_lst = ['ms_th', 'ms2_th', 'hg_th', 'ms_ppm', 'ms2_ppm', 'hg_ppm', 'dda_top', 'sn_ratio',
                          'core_number', 'max_ram', 'img_dpi', 'ms_max']
        f_type_key_lst = ['rt_start', 'rt_end', 'mz_start', 'mz_end', 'pr_window', 'ms2_infopeak_threshold',
                          'ms2_hginfopeak_threshold', 'score_filter', 'isotope_score_filter', 'rank_score_filter']
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
                except configparser.Error as _cfg_e:
                    print('ERROR', _cfg_e)
                    cfg_error += str(_cfg_e)
                    cfg_error += '\n'
        return cfg_params_dct, cfg_error

    def b_run_batchmode(self):

        self.ui.tab_b_statusrun_pte.clear()

        loaded_cfg_files = str(self.ui.tab_b_infiles_pte.toPlainText())
        pre_loaded_cfg_lst = loaded_cfg_files.split('\n')

        max_process = self.ui.tab_b_maxbatch_spb.value()
        sub_max_core = self.ui.tab_b_maxsubcore_spb.value()
        sub_max_ram = self.ui.tab_b_maxsubram_spb.value()

        loaded_cfg_lst = []
        for f in pre_loaded_cfg_lst:
            if len(f) > 4:
                loaded_cfg_lst.append(f)

        tot_num = len(loaded_cfg_lst)
        run_counter = 1

        os.chdir(self.lipidhunter_cwd)

        for _cfg in loaded_cfg_lst:

            self.ui.tab_b_statusrun_pte.insertPlainText('Start processing...\n%s\n' % _cfg)
            hunter_param_dct, cfg_error = self.b_read_cfg(_cfg)
            if len(cfg_error) > 0:
                self.ui.tab_b_statusrun_pte.insertPlainText(str(cfg_error))
            if 'vendor' in list(hunter_param_dct.keys()):
                hunter_param_dct['batch_cfg_file'] = _cfg
                hunter_param_dct['core_number'] = sub_max_core
                hunter_param_dct['max_ram'] = sub_max_ram
                start_time_str = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
                hunter_param_dct['hunter_start_time'] = start_time_str
                os.chdir(hunter_param_dct['hunter_folder'])
                try:
                    os.chdir(hunter_param_dct['hunter_folder'])
                except IOError:
                    print('LipidHunter folder path in configuration is not correct')
                if not os.path.exists(hunter_param_dct['img_output_folder_str']):
                    os.makedirs(hunter_param_dct['img_output_folder_str'])
                param_log_output_path_str = (hunter_param_dct['img_output_folder_str'] +
                                             '/LipidHunter_Params-Log_%s.txt' % hunter_param_dct[
                                                 'hunter_start_time'])
                try:
                    config = configparser.ConfigParser()
                    with open(param_log_output_path_str, 'w') as usr_param_cfg:
                        config.add_section('parameters')
                        for param in list(hunter_param_dct.keys()):
                            config.set('parameters', str(param), str(hunter_param_dct[param]))
                        config.write(usr_param_cfg)
                    log_lst = []
                    hunter_time, log_lst, export_df = huntlipids(hunter_param_dct, error_lst=log_lst)
                    run_time = str(hunter_time)
                    if isinstance(run_time, str):
                        self.ui.tab_b_statusrun_pte.appendPlainText('>>> %s' % run_time)
                        self.ui.tab_b_statusrun_pte.appendPlainText('FINISHED with file %i / %i\n' %
                                                                    (run_counter, tot_num))
                        run_counter += 1
                    else:
                        self.ui.tab_b_statusrun_pte.insertPlainText(
                            '!! Failed to process batch mode configure file:\n Please check settings!!')
                        if len(log_lst) > 0:
                            for err in log_lst:
                                self.ui.tab_b_statusrun_pte.appendPlainText(str(err) + '\n')
                except IOError:
                    self.ui.tab_b_statusrun_pte.appendPlainText('!! Failed to save parameter log files !!')
            else:
                self.ui.tab_b_statusrun_pte.insertPlainText(
                    '!! Failed read batch mode configure files:\n %s \n Please check settings!!' % _cfg)

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
        file_info_str = 'FA white list files (*.xlsx *.XLSX)'
        self.open_file(file_info_str, self.ui.tab_c_fawhitelist_le)

    def c_load_hgcfg(self):
        file_info_str = 'MS Excel files (*.xlsx *.XLSX)'
        self.open_file(file_info_str, self.ui.tab_c_hgcfg_le)

    def c_load_scorecfg(self):
        file_info_str = 'MS Excel files (*.xlsx *.XLSX)'
        self.open_file(file_info_str, self.ui.tab_c_scorecfg_le)


if __name__ == '__main__':
    import sys

    app = QtGui.QApplication(sys.argv)
    window = LipidHunterMain()
    window.show()
    sys.exit(app.exec_())
