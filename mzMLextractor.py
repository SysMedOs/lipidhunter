# -*- coding: utf-8 -*-
# Copyright 2015-2016 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

import os
import glob
import re

from PySide import QtCore, QtGui
import pandas as pd

from mzMLextractor_UI import Ui_MainWindow
from mzMLextractorLib import Extractor
from mzMLextractorLib.Linker import hunt_link


class MainWindow(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None):
        QtGui.QMainWindow.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        # links
        self.ui.logo_lb.setOpenExternalLinks(True)
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
        QtCore.QObject.connect(self.ui.tab_b_addmzml_pb, QtCore.SIGNAL("clicked()"), self.b_load_mzml)
        QtCore.QObject.connect(self.ui.tab_b_addmzmlfolder_pb, QtCore.SIGNAL("clicked()"), self.b_load_mzmlfolder)
        QtCore.QObject.connect(self.ui.tab_b_clearall_pb, QtCore.SIGNAL("clicked()"), self.ui.tab_b_infiles_pte.clear)
        QtCore.QObject.connect(self.ui.tab_b_savexlsxfolder_pb, QtCore.SIGNAL("clicked()"), self.b_save_xls2folder)
        QtCore.QObject.connect(self.ui.tab_b_runextract_pb, QtCore.SIGNAL("clicked()"), self.b_run_extractor)

        # slots for tab d
        QtCore.QObject.connect(self.ui.tab_d_lipidstable_pb, QtCore.SIGNAL("clicked()"), self.d_load_lipidstable)
        QtCore.QObject.connect(self.ui.tab_d_ms2info_pb, QtCore.SIGNAL("clicked()"), self.d_load_ms2info)
        QtCore.QObject.connect(self.ui.tab_d_ms2mzml_pb, QtCore.SIGNAL("clicked()"), self.d_load_mzml)
        QtCore.QObject.connect(self.ui.tab_b_xlsxpath_pb, QtCore.SIGNAL("clicked()"), self.d_save_output)
        QtCore.QObject.connect(self.ui.tab_d_runextract_pb, QtCore.SIGNAL("clicked()"), self.d_run_linker)

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
            print('_mzml', _mzml)
            print('MS_TH', a_ms_th, type(a_ms_th))
            if os.path.isfile(_mzml):
                _mzml_path, _mzml_name = os.path.split(_mzml)
                self.ui.tab_a_statusextractor_pte.insertPlainText(unicode('Start processing...\n%s \n' % _mzml))

                _xlsx_path = _save_xlsx_folder_str + '\\' + _mzml_name[:-4] + 'xlsx'
                _ms_df = extractor.get_ms_all(_mzml, a_ms_th)
                # _ms_df = _ms_df.drop_duplicates(subset=['mz'], keep='first')
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

        cm_pkl_df = cm_pkl_df[['mz', 'i']]
        cm_pkl_df = cm_pkl_df.sort_values(by=['mz', 'i'], ascending=[True, False])
        cm_pkl_df = cm_pkl_df.drop_duplicates(subset=['mz'], keep='first')
        cm_pkl_df = cm_pkl_df.reset_index(drop=True)

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

    def b_load_mzml(self):
        # check existed files
        _loaded_files = str(self.ui.tab_b_infiles_pte.toPlainText())
        _loaded_lst = _loaded_files.split('\n')

        b_load_mzml_dialog = QtGui.QFileDialog(self)
        b_load_mzml_dialog.setNameFilters([u'mzML spectra files (*.mzML *.mzml)'])
        b_load_mzml_dialog.selectNameFilter(u'mzML spectra files (*.mzML *.mzml)')
        if b_load_mzml_dialog.exec_():
            b_load_mzml_str = b_load_mzml_dialog.selectedFiles()[0]
            b_load_mzml_str = os.path.abspath(b_load_mzml_str)
            if b_load_mzml_str not in _loaded_lst:
                self.ui.tab_b_infiles_pte.insertPlainText(unicode(b_load_mzml_str))  # take unicode only
                self.ui.tab_b_infiles_pte.insertPlainText(u'\n')
            else:
                _msgBox = QtGui.QMessageBox()
                _msgBox.setText(u'Spectrum has been chosen already.')
                _msgBox.exec_()

    def b_load_mzmlfolder(self):
        # check existed files
        _loaded_files = str(self.ui.tab_b_infiles_pte.toPlainText())
        _loaded_lst = _loaded_files.split('\n')

        b_load_mzmlfolder_str = QtGui.QFileDialog.getExistingDirectory()
        _mzml_name_lst, _mzml_path_lst = self.get_same_files(b_load_mzmlfolder_str, filetype_lst=['*.mzml', '*.mzML'])
        _duplicated_str = ''
        for _mzml in _mzml_path_lst:
            if _mzml not in _loaded_lst:
                self.ui.tab_b_infiles_pte.insertPlainText(unicode(_mzml))
                self.ui.tab_b_infiles_pte.insertPlainText(u'\n')
            else:
                _duplicated_str = _duplicated_str + unicode(_mzml) + u'\n'
        if len(_duplicated_str) > 0:
            _msgBox = QtGui.QMessageBox()
            _msgBox.setText(_duplicated_str + u'Already chosen. \n Skipped')
            _msgBox.exec_()

    def b_save_xls2folder(self):
        b_save_xlsfolder_str = QtGui.QFileDialog.getExistingDirectory()
        self.ui.tab_b_outpufolder_le.setText(unicode(b_save_xlsfolder_str))
        
    def b_run_extractor(self):
        self.ui.tab_b_statusrun_pte.clear()
        b_ms_th = self.ui.tab_b_msthreshold_spb.value()
        b_ms2_th = self.ui.tab_b_ms2threshold_spb.value()
        self.ui.tab_b_statusrun_pte.insertPlainText(unicode('MS threshold (absolute): %i \n' % b_ms_th))
        extractor = Extractor.Extractor()
        _loaded_mzml_files = str(self.ui.tab_b_infiles_pte.toPlainText())
        _loaded_mzml_lst = _loaded_mzml_files.split('\n')

        _save_xlsx_folder_str = str(self.ui.tab_b_outpufolder_le.text())

        for _mzml in _loaded_mzml_lst:
            if os.path.isfile(_mzml):
                _mzml_path, _mzml_name = os.path.split(_mzml)
                self.ui.tab_b_statusrun_pte.insertPlainText(unicode('Start processing...\n%s \n' % _mzml))

                _xlsx_path = _save_xlsx_folder_str + '\\' + _mzml_name[:-5] + '_all_scan_info.xlsx'
                _xlsx_ms2_path = _save_xlsx_folder_str + '\\' + _mzml_name[:-5] + '_ms2_info.xlsx'
                _ms_df = extractor.get_scan_events(_mzml, b_ms_th, b_ms2_th)
                # _ms_df = _ms_df.drop_duplicates(subset=['mz'], keep='first')
                _ms_df.to_excel(_xlsx_path)
                _ms_df['function'] = _ms_df['function'].apply(pd.to_numeric)
                _ms2_df = _ms_df[_ms_df['function'] > 1]
                _ms2_df = _ms2_df.reset_index()
                _ms2_df.to_excel(_xlsx_ms2_path)
                self.ui.tab_b_statusrun_pte.insertPlainText(unicode('Save as: \n%s.xlsx \n' % _mzml[0:-4]))
        self.ui.tab_b_statusrun_pte.insertPlainText(u'Finished!')

    def d_load_lipidstable(self):
        d_load_lipidstable_dialog = QtGui.QFileDialog(self)
        d_load_lipidstable_dialog.setNameFilters([u'MS Excel files (*.xlsx *.XLSX)'])
        d_load_lipidstable_dialog.selectNameFilter(u'MS Excel files (*.xlsx *.XLSX)')
        if d_load_lipidstable_dialog.exec_():
            self.ui.tab_d_lipidstable_le.clear()
            d_load_mzml_str = d_load_lipidstable_dialog.selectedFiles()[0]
            d_load_mzml_str = os.path.abspath(d_load_mzml_str)
            self.ui.tab_d_lipidstable_le.setText(unicode(d_load_mzml_str))

    def d_load_ms2info(self):
        d_load_ms2info_dialog = QtGui.QFileDialog(self)
        d_load_ms2info_dialog.setNameFilters([u'MS Excel files (*.xlsx *.XLSX)'])
        d_load_ms2info_dialog.selectNameFilter(u'MS Excel files (*.xlsx *.XLSX)')
        if d_load_ms2info_dialog.exec_():
            self.ui.tab_d_ms2info_le.clear()
            d_load_mzml_str = d_load_ms2info_dialog.selectedFiles()[0]
            d_load_mzml_str = os.path.abspath(d_load_mzml_str)
            self.ui.tab_d_ms2info_le.setText(unicode(d_load_mzml_str))

    def d_load_mzml(self):
        d_load_mzml_dialog = QtGui.QFileDialog(self)
        d_load_mzml_dialog.setNameFilters([u'mzML spectra files (*.mzML *.mzml)'])
        d_load_mzml_dialog.selectNameFilter(u'mzML spectra files (*.mzML *.mzml)')
        if d_load_mzml_dialog.exec_():
            self.ui.tab_d_ms2mzml_le.clear()
            d_load_mzml_str = d_load_mzml_dialog.selectedFiles()[0]
            d_load_mzml_str = os.path.abspath(d_load_mzml_str)
            self.ui.tab_d_ms2mzml_le.setText(unicode(d_load_mzml_str))

    def d_save_output(self):
        d_save_output_path = QtGui.QFileDialog.getSaveFileName(caption=u'Save file', filter=u'.xlsx')
        self.ui.tab_d_xlsxpath_le.clear()
        d_load_mzml_str = os.path.abspath(d_save_output_path[0])
        self.ui.tab_d_xlsxpath_le.setText(unicode(d_load_mzml_str))

    def d_run_linker(self):

        print('linker started!')
        _pl_class_info = str(self.ui.tab_d_lipidclass_cmb.currentText())

        pl_class_checker = re.compile(r'(.*)( [\(])(\w{2,3})([\)] )(.*)')

        pl_class_match = pl_class_checker.match(_pl_class_info)

        if pl_class_match:
            pl_class_info_lst = pl_class_match.groups()
            _pl_class = pl_class_info_lst[2]
            _pl_charge = pl_class_info_lst[4]
        else:
            _pl_class = 'PC'
            _pl_charge = '[M+HCOO]-'

        _lipidstable_path_str = str(self.ui.tab_d_lipidstable_le.text())
        _ms2info_path_str = str(self.ui.tab_d_ms2info_le.text())
        _mzml_path_str = str(self.ui.tab_d_ms2mzml_le.text())
        _output_path_str = str(self.ui.tab_d_xlsxpath_le.text())

        ms2_delta = 0.9

        ident_df = pd.read_excel(_lipidstable_path_str, sheetname=0, header=0)
        ms2_df = pd.read_excel(_ms2info_path_str, sheetname=0)

        # obs_mz_lst = ident_df['Input Mass'].tolist()
        ident_idx_lst = ident_df.index.tolist()
        obs_idx_lst = ms2_df.index.tolist()

        step1_df = pd.DataFrame()
        self.ui.tab_d_statusrun_pte.clear()
        self.ui.tab_d_statusrun_pte.insertPlainText(u'Start!')

        for _idx, _ident_se in ident_df.iterrows():

            _obs_mz = _ident_se['Input Mass']
            _lib_mz = _ident_se['Matched Mass']
            _abbr = _ident_se['Abbreviation']
            _formula = _ident_se['Formula']
            _ion = _ident_se['Ion']

            _temp_df = pd.DataFrame()

            _obs_mz_l = _obs_mz - ms2_delta
            _obs_mz_h = _obs_mz + ms2_delta

            _query_code = '%f <= mz <= %f' % (_obs_mz_l, _obs_mz_h)

            _temp_df = ms2_df.query(_query_code)
            if _temp_df.shape[0] > 0:

                _temp_df['MS1_obs_mz'] = _obs_mz
                _temp_df['Lib_mz'] = _lib_mz
                _temp_df['Abbreviation'] = _abbr
                _temp_df['Formula'] = _formula
                _temp_df['Ion'] = _ion
                # print _temp_df
                step1_df = step1_df.append(_temp_df)
            else:
                print _obs_mz, 'not found!'

        _ms_th = self.ui.tab_d_msthreshold_spb.value()
        _ms2_th = self.ui.tab_d_ms2threshold_spb.value()
        _dda_top = self.ui.tab_d_dda_spb.value()
        _rt_start = self.ui.tab_d_rtstart_dspb.value()
        _rt_end = self.ui.tab_d_rtend_dspb.value()
        _mz_start = self.ui.tab_d_mzstart_dspb.value()
        _mz_end = self.ui.tab_d_mzend_dspb.value()

        # Construct parameters for the hunt_link function
        link_params_dct = {'MS_THRESHOLD': _ms_th, 'MS2_THRESHOLD': _ms2_th, 'DDA_TOP': _dda_top,
                           'RT_START': _rt_start, 'RT_END': _rt_end, 'MZ_START': _mz_start, 'MZ_END': _mz_end}

        final_output_df = hunt_link(pl_class=_pl_class, usr_mzml=_mzml_path_str, usr_df=step1_df,
                                    params_dct=link_params_dct)

        final_output_df = final_output_df[final_output_df['MS1_obs_mz'] > 0]

        final_output_df.to_excel(_output_path_str)
        self.ui.tab_d_statusrun_pte.insertPlainText(unicode('Finished!'))


if __name__ == '__main__':
    import sys

    app = QtGui.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
