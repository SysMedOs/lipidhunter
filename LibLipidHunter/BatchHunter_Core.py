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

from __future__ import print_function
from __future__ import division
import ConfigParser as configparser
import glob
import os
import time
import multiprocessing
from multiprocessing import Pool

from PySide import QtCore, QtGui

from BatchHunter_UI import Ui_MainWindow
from LibLipidHunter.SpectraHunter import huntlipids


class BatchHunterCore(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None, cwd=None):
        QtGui.QMainWindow.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        # current folder:
        if cwd is not None:
            print('User LipidHunter folder', cwd)
            self.lipidhunter_cwd = cwd
        else:
            auto_cwd = os.getcwd()
            print('User LipidHunter folder', auto_cwd)
            self.lipidhunter_cwd = auto_cwd

        # links
        self.ui.logo_lb.setOpenExternalLinks(True)

        # slots
        QtCore.QObject.connect(self.ui.tab_a_addcfg_pb, QtCore.SIGNAL("clicked()"), self.a_load_batchcfg)
        QtCore.QObject.connect(self.ui.tab_a_addcfgfolder_pb, QtCore.SIGNAL("clicked()"), self.a_load_batchcfgfolder)
        QtCore.QObject.connect(self.ui.tab_a_clearall_pb, QtCore.SIGNAL("clicked()"), self.ui.tab_b_infiles_pte.clear)
        QtCore.QObject.connect(self.ui.tab_a_runbatch_pb, QtCore.SIGNAL("clicked()"), self.a_run_batch)

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

    @staticmethod
    def read_batch_settings(batch_cfg):
        config = configparser.ConfigParser()
        config.read(batch_cfg)
        batch_cfg_dct = {}
        if config.has_section('parameters'):
            user_cfg = 'parameters'
            options = config.options(user_cfg)
            for param in options:
                batch_cfg_dct[param] = config.get(user_cfg, param)
        batch_cfg_key_lst = batch_cfg_dct.keys()
        i_type_key_lst = ['ms_th', 'ms2_th', 'hg_th', 'ms_ppm', 'ms2_ppm', 'hg_ppm', 'dda_top']
        f_type_key_lst = ['rt_start', 'rt_end', 'mz_start', 'mz_end', 'pr_window',
                          'ms2_infopeak_threshold', 'ms2_hginfopeak_threshold',
                          'score_filter', 'isotope_score_filter']

        if len(batch_cfg_key_lst) > 0:
            for cfg_key in batch_cfg_key_lst:
                if cfg_key in i_type_key_lst:
                    try:
                        batch_cfg_dct[cfg_key] = int(batch_cfg_dct[cfg_key])
                    except ValueError:
                        batch_cfg_dct[cfg_key] = int(float(batch_cfg_dct[cfg_key]))
                elif cfg_key in f_type_key_lst:
                    batch_cfg_dct[cfg_key] = float(batch_cfg_dct[cfg_key])

        return batch_cfg_dct

    def a_load_batchcfg(self):
        # check existed files
        _loaded_files = str(self.ui.tab_b_infiles_pte.toPlainText())
        _loaded_lst = _loaded_files.split('\n')

        b_load_cfg_dialog = QtGui.QFileDialog(self)
        b_load_cfg_dialog.setNameFilters([u'LPPtiger batch mode files (*.txt)'])
        b_load_cfg_dialog.selectNameFilter(u'LPPtiger batch mode files (*.txt)')
        if b_load_cfg_dialog.exec_():
            b_load_cfg_str = b_load_cfg_dialog.selectedFiles()[0]
            b_load_cfg_str = os.path.abspath(b_load_cfg_str)
            if b_load_cfg_str not in _loaded_lst:
                self.ui.tab_b_infiles_pte.insertPlainText(b_load_cfg_str)  # take unicode only
                self.ui.tab_b_infiles_pte.insertPlainText(u'\n')
            else:
                _msgBox = QtGui.QMessageBox()
                _msgBox.setText(u'Batch config file has been chosen already.')
                _msgBox.exec_()

    def a_load_batchcfgfolder(self):
        # check existed files
        _loaded_files = str(self.ui.tab_b_infiles_pte.toPlainText())
        _loaded_lst = _loaded_files.split('\n')

        b_load_cfgfolder_str = QtGui.QFileDialog.getExistingDirectory()
        _cfg_name_lst, _cfg_path_lst = self.get_same_files(b_load_cfgfolder_str, filetype_lst=['*.txt', '*.txt'])
        _duplicated_str = ''
        for _cfg in _cfg_path_lst:
            if _cfg not in _loaded_lst:
                self.ui.tab_b_infiles_pte.insertPlainText(_cfg)
                self.ui.tab_b_infiles_pte.insertPlainText('\n')
            else:
                _duplicated_str = _duplicated_str + _cfg + '\n'
        if len(_duplicated_str) > 0:
            _msgBox = QtGui.QMessageBox()
            _msgBox.setText(_duplicated_str + u'Already chosen. \n Skipped')
            _msgBox.exec_()

    def a_run_batch(self):

        self.ui.tab_b_statusrun_pte.clear()

        loaded_cfg_files = str(self.ui.tab_b_infiles_pte.toPlainText())
        pre_loaded_cfg_lst = loaded_cfg_files.split('\n')

        loaded_cfg_lst = []
        for f in pre_loaded_cfg_lst:
            if len(f) > 4:
                loaded_cfg_lst.append(f)

        tot_num = len(loaded_cfg_lst)
        run_counter = 1

        max_process = self.ui.tab_a_maxbatch_spb.value()
        max_sub_core = self.ui.tab_a_maxsubcore_spb.value()
        max_sub_ram = self.ui.tab_a_maxsubram_spb.value()

        cfg_dct_lst = []
        for cfg_file in loaded_cfg_lst:
            hunter_param_dct = self.read_batch_settings(cfg_file)
            if 'hunter_start_time' in hunter_param_dct.keys():
                hunter_param_dct['batch_cfg_file'] = cfg_file
                cfg_dct_lst.append(hunter_param_dct)
            else:
                hunter_param_dct['batch_cfg_file'] = ''

        if len(cfg_dct_lst) > max_process:
            sub_part_list = map(None, *(iter(cfg_dct_lst),) * max_process)
        else:
            sub_part_list = [cfg_dct_lst]

        tot_part = len(sub_part_list)
        part_num = 1
        for sub_cfg_lst in sub_part_list:
            sub_cfg_lst = filter(lambda x: x is not None, sub_cfg_lst)
            parallel_pool = Pool(max_process)
            hunter_results_lst = []
            core_worker_count = 1
            for _cfg_dct in sub_cfg_lst:
                time.sleep(1)
                start_time_str = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
                _cfg_dct['hunter_start_time'] = start_time_str
                self.ui.tab_b_statusrun_pte.insertPlainText('Start Batch %i / %i file %i / %i ...\n'
                                                            % (part_num, tot_part, core_worker_count, max_process))
                self.ui.tab_b_statusrun_pte.insertPlainText('>>> processing...\n')

                tot_run_time = parallel_pool.apply_async(huntlipids, args=(_cfg_dct,))

                core_worker_count += 1
                hunter_results_lst.append(tot_run_time)

            parallel_pool.close()
            parallel_pool.join()

            for hunter_time in hunter_results_lst:

                run_time = str(hunter_time.get())

                if isinstance(run_time, str):
                    self.ui.tab_b_statusrun_pte.appendPlainText('>>> %s sec ' % run_time)
                    self.ui.tab_b_statusrun_pte.appendPlainText('FINISHED with file %i / %i\n\n' %
                                                                (run_counter, tot_num))
                    run_counter += 1
                else:
                    self.ui.tab_b_statusrun_pte.insertPlainText(
                        '!! Failed to process batch mode configure file:\n Please check settings!!')

            part_num += 1


if __name__ == '__main__':
    multiprocessing.freeze_support()
    import sys
    app = QtGui.QApplication(sys.argv)
    window = BatchHunterCore()
    window.show()
    sys.exit(app.exec_())
