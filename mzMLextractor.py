# -*- coding: utf-8 -*-
# Copyright 2015-2016 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

from __future__ import print_function

import sys
import os
import glob
from PySide import QtCore, QtGui

from mzMLextractor_UI import Ui_MainWindow


class MainWindow(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None):
        QtGui.QMainWindow.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        # slots for tab a
        QtCore.QObject.connect(self.ui.tab_a_addmzml_pb, QtCore.SIGNAL("clicked()"), self.a_load_mzml)
        QtCore.QObject.connect(self.ui.tab_a_addmzmlfolder_pb, QtCore.SIGNAL("clicked()"), self.a_load_mzmlfolder)
        QtCore.QObject.connect(self.ui.tab_a_clearall_pb, QtCore.SIGNAL("clicked()"), self.ui.tab_a_infiles_pte.clear)

        # slots for tab b
        QtCore.QObject.connect(self.ui.tab_b_clearall_pb, QtCore.SIGNAL("clicked()"), self.ui.tab_b_infiles_pte.clear)

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
        if a_load_mzmlfolder_str is not u'':
            os.chdir(a_load_mzmlfolder_str)
            # capability for Linux and Windows
            _mzml_lower_lst = glob.glob('*.mzml')
            _mzml_upper_lst = glob.glob('*.mzML')
            # merge list
            _pre_mzml_lst = _mzml_lower_lst + [m for m in _mzml_upper_lst if m not in _mzml_lower_lst]
            _mzml_lst = list(os.path.abspath(m) for m in _pre_mzml_lst)
            _duplicated_str = ''
            for _mzml in _mzml_lst:
                if _mzml not in _loaded_lst:
                    self.ui.tab_a_infiles_pte.insertPlainText(unicode(_mzml))
                    self.ui.tab_a_infiles_pte.insertPlainText(u'\n')
                else:
                    _duplicated_str = _duplicated_str + unicode(_mzml) + u'\n'
            if len(_duplicated_str) > 0:
                _msgBox = QtGui.QMessageBox()
                _msgBox.setText(_duplicated_str + u'Already chosen. \n Skipped')
                _msgBox.exec_()


if __name__ == '__main__':
    import sys

    app = QtGui.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
