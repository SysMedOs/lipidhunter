# -*- coding: utf-8 -*-
# Copyright 2015-2016 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

import sys
import os
import glob

from PySide import QtCore, QtGui
import pandas as pd

from mzMLextractor_UI import Ui_MainWindow
from mzMLextractorLib import Extractor


class MainWindow(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None):
        QtGui.QMainWindow.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        # links
        self.ui.label_5.setOpenExternalLinks(True)
        self.ui.tab_c_3_lb.setOpenExternalLinks(True)
        self.ui.tab_c_5_lb.setOpenExternalLinks(True)

        # slots for tab a
        QtCore.QObject.connect(self.ui.tab_a_addmzml_pb, QtCore.SIGNAL("clicked()"), self.a_load_mzml)
        QtCore.QObject.connect(self.ui.tab_a_addmzmlfolder_pb, QtCore.SIGNAL("clicked()"), self.a_load_mzmlfolder)
        QtCore.QObject.connect(self.ui.tab_a_clearall_pb, QtCore.SIGNAL("clicked()"), self.ui.tab_a_infiles_pte.clear)
        QtCore.QObject.connect(self.ui.tab_a_savexlsxfolder_pb, QtCore.SIGNAL("clicked()"), self.a_save_xls2folder)
        QtCore.QObject.connect(self.ui.tab_a_savecsv_pb, QtCore.SIGNAL("clicked()"), self.a_save_csv2folder)
        QtCore.QObject.connect(self.ui.tab_a_runextractor_pb, QtCore.SIGNAL("clicked()"), self.a_run_extractor)
        QtCore.QObject.connect(self.ui.tab_a_runmerge_pb, QtCore.SIGNAL("clicked()"), self.a_run_merger)

        # slots for tab b
        QtCore.QObject.connect(self.ui.tab_b_clearall_pb, QtCore.SIGNAL("clicked()"), self.ui.tab_b_infiles_pte.clear)

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
                _pre_found_lst = _pre_found_lst + [f for f in _tmp_found_lst if f not in _pre_found_lst]
            filename_lst = _pre_found_lst
            abs_path_lst = list(os.path.abspath(ff) for ff in _pre_found_lst)
        else:
            filename_lst = []
            abs_path_lst = []

        return filename_lst, abs_path_lst

    def a_load_mzml(self):
        # check existed files
        _loaded_files = str(self.ui.tab_a_infiles_pte.toPlainText())
        _loaded_lst = _loaded_files.split('\n')

        a_load_mzml_dialog = QtGui.QFileDialog(self)
        a_load_mzml_dialog.setNameFilters([u'mzML spectra files (*.mzML *.mzml)'])
        a_load_mzml_dialog.selectNameFilter(u'mzML spectra files (*.mzML *.mzml)')
        if a_load_mzml_dialog.exec_():
            a_load_mzml_str = a_load_mzml_dialog.selectedFiles()[0]
            a_load_mzml_str = os.path.abspath(a_load_mzml_str)
            if a_load_mzml_str not in _loaded_lst:
                self.ui.tab_a_infiles_pte.insertPlainText(unicode(a_load_mzml_str))  # take unicode only
                self.ui.tab_a_infiles_pte.insertPlainText(u'\n')
            else:
                _msgBox = QtGui.QMessageBox()
                _msgBox.setText(u'Spectrum has been chosen already.')
                _msgBox.exec_()

    def a_load_mzmlfolder(self):
        # check existed files
        _loaded_files = str(self.ui.tab_a_infiles_pte.toPlainText())
        _loaded_lst = _loaded_files.split('\n')

        a_load_mzmlfolder_str = QtGui.QFileDialog.getExistingDirectory()
        _mzml_name_lst, _mzml_path_lst = self.get_same_files(a_load_mzmlfolder_str, filetype_lst=['*.mzml', '*.mzML'])
        _duplicated_str = ''
        for _mzml in _mzml_path_lst:
            if _mzml not in _loaded_lst:
                self.ui.tab_a_infiles_pte.insertPlainText(unicode(_mzml))
                self.ui.tab_a_infiles_pte.insertPlainText(u'\n')
            else:
                _duplicated_str = _duplicated_str + unicode(_mzml) + u'\n'
        if len(_duplicated_str) > 0:
            _msgBox = QtGui.QMessageBox()
            _msgBox.setText(_duplicated_str + u'Already chosen. \n Skipped')
            _msgBox.exec_()


    def a_save_xls2folder(self):
        a_save_xlsfolder_str = QtGui.QFileDialog.getExistingDirectory()
        self.ui.tab_a_xlsxfolder_le.setText(unicode(a_save_xlsfolder_str))

    def a_save_csv2folder(self):
        a_save_csvfolder_str = QtGui.QFileDialog.getSaveFileName(caption=u'Save file', filter=u'.csv')
        self.ui.tab_a_csvfolder_le.setText(unicode(a_save_csvfolder_str[0]))

    def a_run_extractor(self):
        self.ui.tab_a_statusextractor_pte.clear()
        a_ms_th = self.ui.tab_a_msthreshold_spb.value()
        self.ui.tab_a_statusextractor_pte.insertPlainText(unicode('MS threshold (absolute): %i \n' % a_ms_th))
        extractor = Extractor.Extractor()
        _loaded_mzml_files = str(self.ui.tab_a_infiles_pte.toPlainText())
        _loaded_mzml_lst = _loaded_mzml_files.split('\n')

        _save_xlsx_folder_str = str(self.ui.tab_a_xlsxfolder_le.text())

        for _mzml in _loaded_mzml_lst:
            if os.path.isfile(_mzml):
                _mzml_path, _mzml_name = os.path.split(_mzml)
                self.ui.tab_a_statusextractor_pte.insertPlainText(unicode('Start processing...\n%s \n' % _mzml))

                _xlsx_path = _save_xlsx_folder_str + '\\' + _mzml_name[:-4] + 'xlsx'
                _ms_df = extractor.get_ms_all(_mzml, a_ms_th)
                _ms_df.to_excel(_xlsx_path)
                self.ui.tab_a_statusextractor_pte.insertPlainText(unicode('Save as: \n%s.xlsx \n' % _mzml[0:-4]))
        self.ui.tab_a_statusextractor_pte.insertPlainText(u'Finished!')

    def a_run_merger(self):
        self.ui.tab_a_statusmerger_pte.insertPlainText(u'Start to proceed...')

        _save_csv_str = str(self.ui.tab_a_csvfolder_le.text())

        _save_xlsx_folder_str = str(self.ui.tab_a_xlsxfolder_le.text())

        _xlsx_name_lst, _xlsx_path_lst = self.get_same_files(_save_xlsx_folder_str, filetype_lst=['*.xlsx', '*.XLSX'])
        cm_pkl_lst = zip(_xlsx_name_lst, _xlsx_path_lst)

        cm_pkl_df = pd.DataFrame()
        for _cm in cm_pkl_lst:
            self.ui.tab_a_statusmerger_pte.insertPlainText(unicode('reading --> %s \n' % _cm[0]))
            _tmp_df = pd.read_excel(_cm[1])
            _tmp_df['file'] = _cm[0]
            cm_pkl_df = cm_pkl_df.append(_tmp_df)

        cm_pkl_df['mz_2f'] = cm_pkl_df['mz']
        cm_pkl_df['rt_2f'] = cm_pkl_df['rt']
        cm_pkl_df = cm_pkl_df.round({'mz_2f': 2, 'rt_2f': 2})
        cm_pkl_df = cm_pkl_df.sort_values(by=['mz', 'rt'])
        cm_pkl_df = cm_pkl_df.reset_index()
        if cm_pkl_df.shape[0] > 500000:
            cm_pkl_df_p1 = cm_pkl_df[:500000, :]
            print cm_pkl_df_p1.shape
            cm_pkl_df_p1.to_csv(''.join([_save_csv_str[0:-4], ['_1.csv']]))
            cm_pkl_df_p2 = cm_pkl_df[500000:, :]
            print cm_pkl_df_p2.shape
            cm_pkl_df_p2.to_csv(''.join([_save_csv_str[0:-4], ['_2.csv']]))
            _save_csv_str = ''.join([_save_csv_str[0:-4], ['_1.csv'], '\n', _save_csv_str[0:-4], ['_2.csv']])
        else:
            cm_pkl_df.to_csv(_save_csv_str)

        self.ui.tab_a_statusmerger_pte.insertPlainText(unicode('Merged and saved as %s' % _save_csv_str))

if __name__ == '__main__':
    import sys

    app = QtGui.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
