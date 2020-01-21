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

# Form implementation generated from reading ui file '.\LipidHunter_UI.ui'
#
# Created: Fri May 11 16:51:42 2018
#      by: pyside-uic 0.2.15 running on PySide 1.2.4
# Modified: Tue Jan 21 2020
#
# WARNING! All changes made in this file will be lost!

from PySide2 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow, scale=1):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(960 * scale, 720 * scale)
        icon = QtGui.QIcon()
        icon.addPixmap(
            QtGui.QPixmap(":/LipidHunter.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off
        )
        MainWindow.setWindowIcon(icon)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.tabframe = QtWidgets.QTabWidget(self.centralwidget)
        self.tabframe.setGeometry(QtCore.QRect(10, 10, 931 * scale, 641 * scale))
        self.tabframe.setObjectName("tabframe")
        self.hunter_tab = QtWidgets.QWidget()
        self.hunter_tab.setObjectName("hunter_tab")
        self.gridLayoutWidget_5 = QtWidgets.QWidget(self.hunter_tab)
        self.gridLayoutWidget_5.setGeometry(
            QtCore.QRect(10, 50 * scale, 901 * scale, 341 * scale)
        )
        self.gridLayoutWidget_5.setObjectName("gridLayoutWidget_5")
        self.gridLayout_5 = QtWidgets.QGridLayout(self.gridLayoutWidget_5)
        self.gridLayout_5.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.tab_a_saveimgfolder_le = QtWidgets.QLineEdit(self.gridLayoutWidget_5)
        self.tab_a_saveimgfolder_le.setObjectName("tab_a_saveimgfolder_le")
        self.gridLayout_5.addWidget(self.tab_a_saveimgfolder_le, 5, 1, 1, 1)
        self.label_50 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.label_50.setObjectName("label_50")
        self.gridLayout_5.addWidget(self.label_50, 2, 0, 1, 1)
        self.tab_a_mzml_le = QtWidgets.QLineEdit(self.gridLayoutWidget_5)
        self.tab_a_mzml_le.setObjectName("tab_a_mzml_le")
        self.gridLayout_5.addWidget(self.tab_a_mzml_le, 3, 1, 1, 1)
        self.line_6 = QtWidgets.QFrame(self.gridLayoutWidget_5)
        self.line_6.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_6.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_6.setObjectName("line_6")
        self.gridLayout_5.addWidget(self.line_6, 7, 1, 1, 1)
        self.tab_a_mzml_pb = QtWidgets.QPushButton(self.gridLayoutWidget_5)
        self.tab_a_mzml_pb.setObjectName("tab_a_mzml_pb")
        self.gridLayout_5.addWidget(self.tab_a_mzml_pb, 3, 2, 1, 1)
        self.tab_a_savexlsxpath_le = QtWidgets.QLineEdit(self.gridLayoutWidget_5)
        self.tab_a_savexlsxpath_le.setObjectName("tab_a_savexlsxpath_le")
        self.gridLayout_5.addWidget(self.tab_a_savexlsxpath_le, 6, 1, 1, 1)
        self.label_3 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.label_3.setObjectName("label_3")
        self.gridLayout_5.addWidget(self.label_3, 0, 0, 1, 1)
        self.horizontalLayout_11 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_11.setObjectName("horizontalLayout_11")
        self.label_4 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.label_4.setObjectName("label_4")
        self.horizontalLayout_11.addWidget(self.label_4)
        self.tab_a_rtstart_dspb = QtWidgets.QDoubleSpinBox(self.gridLayoutWidget_5)
        self.tab_a_rtstart_dspb.setProperty("value", 10.0)
        self.tab_a_rtstart_dspb.setObjectName("tab_a_rtstart_dspb")
        self.horizontalLayout_11.addWidget(self.tab_a_rtstart_dspb)
        self.label_9 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.label_9.setObjectName("label_9")
        self.horizontalLayout_11.addWidget(self.label_9)
        self.tab_a_rtend_dspb = QtWidgets.QDoubleSpinBox(self.gridLayoutWidget_5)
        self.tab_a_rtend_dspb.setProperty("value", 25.0)
        self.tab_a_rtend_dspb.setObjectName("tab_a_rtend_dspb")
        self.horizontalLayout_11.addWidget(self.tab_a_rtend_dspb)
        spacerItem = QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_11.addItem(spacerItem)
        self.label_10 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.label_10.setObjectName("label_10")
        self.horizontalLayout_11.addWidget(self.label_10)
        self.tab_a_mzstart_dspb = QtWidgets.QDoubleSpinBox(self.gridLayoutWidget_5)
        self.tab_a_mzstart_dspb.setDecimals(4)
        self.tab_a_mzstart_dspb.setMaximum(1200.0)
        self.tab_a_mzstart_dspb.setProperty("value", 600.0)
        self.tab_a_mzstart_dspb.setObjectName("tab_a_mzstart_dspb")
        self.horizontalLayout_11.addWidget(self.tab_a_mzstart_dspb)
        self.label_11 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.label_11.setObjectName("label_11")
        self.horizontalLayout_11.addWidget(self.label_11)
        self.tab_a_mzend_dspb = QtWidgets.QDoubleSpinBox(self.gridLayoutWidget_5)
        self.tab_a_mzend_dspb.setDecimals(4)
        self.tab_a_mzend_dspb.setMaximum(2000.0)
        self.tab_a_mzend_dspb.setProperty("value", 1000.0)
        self.tab_a_mzend_dspb.setObjectName("tab_a_mzend_dspb")
        self.horizontalLayout_11.addWidget(self.tab_a_mzend_dspb)
        spacerItem1 = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_11.addItem(spacerItem1)
        self.gridLayout_5.addLayout(self.horizontalLayout_11, 8, 1, 1, 1)
        self.label_6 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.label_6.setObjectName("label_6")
        self.gridLayout_5.addWidget(self.label_6, 3, 0, 1, 1)
        self.tab_a_loadfalist_pb = QtWidgets.QPushButton(self.gridLayoutWidget_5)
        self.tab_a_loadfalist_pb.setObjectName("tab_a_loadfalist_pb")
        self.gridLayout_5.addWidget(self.tab_a_loadfalist_pb, 1, 2, 1, 1)
        self.label_7 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.label_7.setObjectName("label_7")
        self.gridLayout_5.addWidget(self.label_7, 5, 0, 1, 1)
        self.tab_a_loadscorecfg_pb = QtWidgets.QPushButton(self.gridLayoutWidget_5)
        self.tab_a_loadscorecfg_pb.setObjectName("tab_a_loadscorecfg_pb")
        self.gridLayout_5.addWidget(self.tab_a_loadscorecfg_pb, 2, 2, 1, 1)
        self.tab_a_sumxlsxpath_pb = QtWidgets.QPushButton(self.gridLayoutWidget_5)
        self.tab_a_sumxlsxpath_pb.setObjectName("tab_a_sumxlsxpath_pb")
        self.gridLayout_5.addWidget(self.tab_a_sumxlsxpath_pb, 6, 2, 1, 1)
        self.tab_a_loadscorecfg_le = QtWidgets.QLineEdit(self.gridLayoutWidget_5)
        self.tab_a_loadscorecfg_le.setObjectName("tab_a_loadscorecfg_le")
        self.gridLayout_5.addWidget(self.tab_a_loadscorecfg_le, 2, 1, 1, 1)
        self.tab_a_saveimgfolder_pb = QtWidgets.QPushButton(self.gridLayoutWidget_5)
        self.tab_a_saveimgfolder_pb.setObjectName("tab_a_saveimgfolder_pb")
        self.gridLayout_5.addWidget(self.tab_a_saveimgfolder_pb, 5, 2, 1, 1)
        self.label_8 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.label_8.setObjectName("label_8")
        self.gridLayout_5.addWidget(self.label_8, 6, 0, 1, 1)
        self.line_5 = QtWidgets.QFrame(self.gridLayoutWidget_5)
        self.line_5.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_5.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_5.setObjectName("line_5")
        self.gridLayout_5.addWidget(self.line_5, 4, 1, 1, 1)
        self.horizontalLayout_15 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_15.setObjectName("horizontalLayout_15")
        self.label_16 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.label_16.setObjectName("label_16")
        self.horizontalLayout_15.addWidget(self.label_16)
        self.tab_a_ms2ppm_spb = QtWidgets.QSpinBox(self.gridLayoutWidget_5)
        self.tab_a_ms2ppm_spb.setMaximum(9999)
        self.tab_a_ms2ppm_spb.setProperty("value", 50)
        self.tab_a_ms2ppm_spb.setObjectName("tab_a_ms2ppm_spb")
        self.horizontalLayout_15.addWidget(self.tab_a_ms2ppm_spb)
        spacerItem2 = QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_15.addItem(spacerItem2)
        self.tab_d_ms2threshold_spb_2 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.tab_d_ms2threshold_spb_2.setObjectName("tab_d_ms2threshold_spb_2")
        self.horizontalLayout_15.addWidget(self.tab_d_ms2threshold_spb_2)
        self.tab_a_ms2threshold_spb = QtWidgets.QSpinBox(self.gridLayoutWidget_5)
        self.tab_a_ms2threshold_spb.setMaximum(999999)
        self.tab_a_ms2threshold_spb.setProperty("value", 10)
        self.tab_a_ms2threshold_spb.setObjectName("tab_a_ms2threshold_spb")
        self.horizontalLayout_15.addWidget(self.tab_a_ms2threshold_spb)
        spacerItem3 = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_15.addItem(spacerItem3)
        self.gridLayout_5.addLayout(self.horizontalLayout_15, 11, 1, 1, 1)
        self.horizontalLayout_7 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_7.setObjectName("horizontalLayout_7")
        self.tab_a_lipidclass_cmb = QtWidgets.QComboBox(self.gridLayoutWidget_5)
        self.tab_a_lipidclass_cmb.setMinimumSize(QtCore.QSize(220, 0))
        self.tab_a_lipidclass_cmb.setMaxVisibleItems(14)
        self.tab_a_lipidclass_cmb.setObjectName("tab_a_lipidclass_cmb")
        self.tab_a_lipidclass_cmb.addItem("")
        self.tab_a_lipidclass_cmb.addItem("")
        self.tab_a_lipidclass_cmb.addItem("")
        self.tab_a_lipidclass_cmb.addItem("")
        self.tab_a_lipidclass_cmb.addItem("")
        self.tab_a_lipidclass_cmb.addItem("")
        self.tab_a_lipidclass_cmb.addItem("")
        self.tab_a_lipidclass_cmb.addItem("")
        self.tab_a_lipidclass_cmb.addItem("")
        self.tab_a_lipidclass_cmb.addItem("")
        self.tab_a_lipidclass_cmb.addItem("")
        self.tab_a_lipidclass_cmb.addItem("")
        self.tab_a_lipidclass_cmb.addItem("")
        self.tab_a_lipidclass_cmb.addItem("")
        self.tab_a_lipidclass_cmb.addItem("")
        self.tab_a_lipidclass_cmb.addItem("")
        self.tab_a_lipidclass_cmb.addItem("")
        self.tab_a_lipidclass_cmb.addItem("")
        self.tab_a_lipidclass_cmb.addItem("")
        # self.tab_a_lipidclass_cmb.addItem("")
        self.horizontalLayout_7.addWidget(self.tab_a_lipidclass_cmb)
        spacerItem4 = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_7.addItem(spacerItem4)
        self.gridLayout_5.addLayout(self.horizontalLayout_7, 0, 1, 1, 1)
        self.label_5 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.label_5.setObjectName("label_5")
        self.gridLayout_5.addWidget(self.label_5, 1, 0, 1, 1)
        self.horizontalLayout_28 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_28.setObjectName("horizontalLayout_28")
        self.tab_a_loadfalist_le = QtWidgets.QLineEdit(self.gridLayoutWidget_5)
        self.tab_a_loadfalist_le.setObjectName("tab_a_loadfalist_le")
        self.horizontalLayout_28.addWidget(self.tab_a_loadfalist_le)
        self.gridLayout_5.addLayout(self.horizontalLayout_28, 1, 1, 1, 1)
        self.horizontalLayout_13 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_13.setObjectName("horizontalLayout_13")
        self.label_19 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.label_19.setObjectName("label_19")
        self.horizontalLayout_13.addWidget(self.label_19)
        self.tab_a_isotopescore_spb = QtWidgets.QDoubleSpinBox(self.gridLayoutWidget_5)
        self.tab_a_isotopescore_spb.setMinimumSize(QtCore.QSize(0, 0))
        self.tab_a_isotopescore_spb.setPrefix("")
        self.tab_a_isotopescore_spb.setDecimals(1)
        self.tab_a_isotopescore_spb.setProperty("value", 80.0)
        self.tab_a_isotopescore_spb.setObjectName("tab_a_isotopescore_spb")
        self.horizontalLayout_13.addWidget(self.tab_a_isotopescore_spb)
        spacerItem5 = QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_13.addItem(spacerItem5)
        self.label_20 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.label_20.setObjectName("label_20")
        self.horizontalLayout_13.addWidget(self.label_20)
        self.tab_a_score_spb = QtWidgets.QDoubleSpinBox(self.gridLayoutWidget_5)
        self.tab_a_score_spb.setMinimumSize(QtCore.QSize(0, 0))
        self.tab_a_score_spb.setPrefix("")
        self.tab_a_score_spb.setDecimals(1)
        self.tab_a_score_spb.setProperty("value", 40.0)
        self.tab_a_score_spb.setObjectName("tab_a_score_spb")
        self.horizontalLayout_13.addWidget(self.tab_a_score_spb)
        spacerItem6 = QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_13.addItem(spacerItem6)
        self.label_38 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.label_38.setObjectName("label_38")
        self.horizontalLayout_13.addWidget(self.label_38)
        self.tab_a_ms2infoth_dspb = QtWidgets.QDoubleSpinBox(self.gridLayoutWidget_5)
        self.tab_a_ms2infoth_dspb.setMaximumSize(QtCore.QSize(75, 16777215))
        self.tab_a_ms2infoth_dspb.setProperty("value", 0.1)
        self.tab_a_ms2infoth_dspb.setObjectName("tab_a_ms2infoth_dspb")
        self.horizontalLayout_13.addWidget(self.tab_a_ms2infoth_dspb)
        spacerItem7 = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_13.addItem(spacerItem7)
        self.gridLayout_5.addLayout(self.horizontalLayout_13, 12, 1, 1, 1)
        self.horizontalLayout_12 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_12.setObjectName("horizontalLayout_12")
        self.label_15 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.label_15.setObjectName("label_15")
        self.horizontalLayout_12.addWidget(self.label_15)
        self.tab_a_msppm_spb = QtWidgets.QSpinBox(self.gridLayoutWidget_5)
        self.tab_a_msppm_spb.setMaximum(9999)
        self.tab_a_msppm_spb.setProperty("value", 20)
        self.tab_a_msppm_spb.setObjectName("tab_a_msppm_spb")
        self.horizontalLayout_12.addWidget(self.tab_a_msppm_spb)
        spacerItem8 = QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_12.addItem(spacerItem8)
        self.label_13 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.label_13.setObjectName("label_13")
        self.horizontalLayout_12.addWidget(self.label_13)
        self.tab_a_msthreshold_spb = QtWidgets.QSpinBox(self.gridLayoutWidget_5)
        self.tab_a_msthreshold_spb.setMaximum(999999)
        self.tab_a_msthreshold_spb.setProperty("value", 1000)
        self.tab_a_msthreshold_spb.setObjectName("tab_a_msthreshold_spb")
        self.horizontalLayout_12.addWidget(self.tab_a_msthreshold_spb)
        spacerItem9 = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_12.addItem(spacerItem9)
        self.tab_a_msmax_chb = QtWidgets.QCheckBox(self.gridLayoutWidget_5)
        self.tab_a_msmax_chb.setObjectName("tab_a_msmax_chb")
        self.horizontalLayout_12.addWidget(self.tab_a_msmax_chb)
        self.gridLayout_5.addLayout(self.horizontalLayout_12, 10, 1, 1, 1)
        self.horizontalLayout_16 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_16.setObjectName("horizontalLayout_16")
        self.label_25 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.label_25.setObjectName("label_25")
        self.horizontalLayout_16.addWidget(self.label_25)
        self.tab_a_prwindow_spb = QtWidgets.QDoubleSpinBox(self.gridLayoutWidget_5)
        self.tab_a_prwindow_spb.setMinimumSize(QtCore.QSize(70, 0))
        self.tab_a_prwindow_spb.setMaximum(2.0)
        self.tab_a_prwindow_spb.setSingleStep(0.01)
        self.tab_a_prwindow_spb.setProperty("value", 0.75)
        self.tab_a_prwindow_spb.setObjectName("tab_a_prwindow_spb")
        self.horizontalLayout_16.addWidget(self.tab_a_prwindow_spb)
        spacerItem10 = QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_16.addItem(spacerItem10)
        self.x = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.x.setObjectName("x")
        self.horizontalLayout_16.addWidget(self.x)
        self.tab_a_dda_spb = QtWidgets.QSpinBox(self.gridLayoutWidget_5)
        self.tab_a_dda_spb.setMaximum(20)
        self.tab_a_dda_spb.setProperty("value", 6)
        self.tab_a_dda_spb.setObjectName("tab_a_dda_spb")
        self.horizontalLayout_16.addWidget(self.tab_a_dda_spb)
        spacerItem11 = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_16.addItem(spacerItem11)
        self.gridLayout_5.addLayout(self.horizontalLayout_16, 9, 1, 1, 1)
        self.tab_a_msmax_spb = QtWidgets.QSpinBox(self.gridLayoutWidget_5)
        self.tab_a_msmax_spb.setMaximum(999999)
        self.tab_a_msmax_spb.setProperty("value", 1000)
        self.tab_a_msmax_spb.setObjectName("tab_a_msmax_spb")
        self.gridLayout_5.addWidget(self.tab_a_msmax_spb, 10, 2, 1, 1)
        spacerItem12 = QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
        )
        self.gridLayout_5.addItem(spacerItem12, 8, 0, 1, 1)
        spacerItem13 = QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
        )
        self.gridLayout_5.addItem(spacerItem13, 10, 0, 1, 1)
        spacerItem14 = QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
        )
        self.gridLayout_5.addItem(spacerItem14, 9, 0, 1, 1)
        spacerItem15 = QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
        )
        self.gridLayout_5.addItem(spacerItem15, 11, 0, 1, 1)
        spacerItem16 = QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
        )
        self.gridLayout_5.addItem(spacerItem16, 12, 0, 1, 1)
        self.runhunter_tabframe = QtWidgets.QTabWidget(self.hunter_tab)
        self.runhunter_tabframe.setGeometry(
            QtCore.QRect(10, 420 * scale, 901 * scale, 181 * scale)
        )
        self.runhunter_tabframe.setObjectName("runhunter_tabframe")
        self.singlerun_tab = QtWidgets.QWidget()
        self.singlerun_tab.setObjectName("singlerun_tab")
        self.horizontalLayoutWidget_3 = QtWidgets.QWidget(self.singlerun_tab)
        self.horizontalLayoutWidget_3.setGeometry(
            QtCore.QRect(169 * scale, 9, 721 * scale, 131 * scale)
        )
        self.horizontalLayoutWidget_3.setObjectName("horizontalLayoutWidget_3")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_3)
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.tab_a_statusrun_pte = QtWidgets.QPlainTextEdit(
            self.horizontalLayoutWidget_3
        )
        self.tab_a_statusrun_pte.setObjectName("tab_a_statusrun_pte")
        self.horizontalLayout.addWidget(self.tab_a_statusrun_pte)
        self.verticalLayoutWidget_8 = QtWidgets.QWidget(self.singlerun_tab)
        self.verticalLayoutWidget_8.setGeometry(
            QtCore.QRect(10, 10, 152 * scale, 129 * scale)
        )
        self.verticalLayoutWidget_8.setObjectName("verticalLayoutWidget_8")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_8)
        self.verticalLayout_5.setSizeConstraint(QtWidgets.QLayout.SetNoConstraint)
        self.verticalLayout_5.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.tab_a_runhunter_pb = QtWidgets.QPushButton(self.verticalLayoutWidget_8)
        self.tab_a_runhunter_pb.setMinimumSize(QtCore.QSize(150, 0))
        self.tab_a_runhunter_pb.setMaximumSize(QtCore.QSize(170, 16777215))
        self.tab_a_runhunter_pb.setObjectName("tab_a_runhunter_pb")
        self.verticalLayout_5.addWidget(self.tab_a_runhunter_pb)
        self.tab_a_runhunter_pgb = QtWidgets.QProgressBar(self.verticalLayoutWidget_8)
        self.tab_a_runhunter_pgb.setMaximumSize(QtCore.QSize(150, 16777215))
        self.tab_a_runhunter_pgb.setProperty("value", 5)
        self.tab_a_runhunter_pgb.setTextVisible(False)
        self.tab_a_runhunter_pgb.setObjectName("tab_a_runhunter_pgb")
        self.verticalLayout_5.addWidget(self.tab_a_runhunter_pgb)
        self.runhunter_tabframe.addTab(self.singlerun_tab, "")
        self.savecfg_tab = QtWidgets.QWidget()
        self.savecfg_tab.setObjectName("savecfg_tab")
        self.verticalLayoutWidget_2 = QtWidgets.QWidget(self.savecfg_tab)
        self.verticalLayoutWidget_2.setGeometry(
            QtCore.QRect(10, 5, 881 * scale, 131 * scale)
        )
        self.verticalLayoutWidget_2.setObjectName("verticalLayoutWidget_2")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_2)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label = QtWidgets.QLabel(self.verticalLayoutWidget_2)
        self.label.setObjectName("label")
        self.horizontalLayout_2.addWidget(self.label)
        self.tab_a_cfgpath_le = QtWidgets.QLineEdit(self.verticalLayoutWidget_2)
        self.tab_a_cfgpath_le.setObjectName("tab_a_cfgpath_le")
        self.horizontalLayout_2.addWidget(self.tab_a_cfgpath_le)
        self.tab_a_cfgpath_pb = QtWidgets.QPushButton(self.verticalLayoutWidget_2)
        self.tab_a_cfgpath_pb.setObjectName("tab_a_cfgpath_pb")
        self.horizontalLayout_2.addWidget(self.tab_a_cfgpath_pb)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.tab_a_gencfg_pb = QtWidgets.QPushButton(self.verticalLayoutWidget_2)
        self.tab_a_gencfg_pb.setObjectName("tab_a_gencfg_pb")
        self.horizontalLayout_3.addWidget(self.tab_a_gencfg_pb)
        self.tab_a_gencfg_pte = QtWidgets.QPlainTextEdit(self.verticalLayoutWidget_2)
        self.tab_a_gencfg_pte.setObjectName("tab_a_gencfg_pte")
        self.horizontalLayout_3.addWidget(self.tab_a_gencfg_pte)
        self.verticalLayout.addLayout(self.horizontalLayout_3)
        self.runhunter_tabframe.addTab(self.savecfg_tab, "")
        self.horizontalLayoutWidget_19 = QtWidgets.QWidget(self.hunter_tab)
        self.horizontalLayoutWidget_19.setGeometry(
            QtCore.QRect(10, 10, 901 * scale, 31 * scale)
        )
        self.horizontalLayoutWidget_19.setObjectName("horizontalLayoutWidget_19")
        self.horizontalLayout_31 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_19)
        self.horizontalLayout_31.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_31.setObjectName("horizontalLayout_31")
        self.horizontalLayout_14 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_14.setObjectName("horizontalLayout_14")
        self.label_23 = QtWidgets.QLabel(self.horizontalLayoutWidget_19)
        self.label_23.setMaximumSize(QtCore.QSize(120, 16777215))
        self.label_23.setObjectName("label_23")
        self.horizontalLayout_14.addWidget(self.label_23)
        self.vendor_cmb = QtWidgets.QComboBox(self.horizontalLayoutWidget_19)
        self.vendor_cmb.setMinimumSize(QtCore.QSize(250, 0))
        self.vendor_cmb.setMaxVisibleItems(15)
        self.vendor_cmb.setObjectName("vendor_cmb")
        self.vendor_cmb.addItem("")
        self.vendor_cmb.addItem("")
        self.vendor_cmb.addItem("")
        self.vendor_cmb.addItem("")
        self.vendor_cmb.addItem("")
        self.horizontalLayout_14.addWidget(self.vendor_cmb)
        self.line_9 = QtWidgets.QFrame(self.horizontalLayoutWidget_19)
        self.line_9.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_9.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_9.setObjectName("line_9")
        self.horizontalLayout_14.addWidget(self.line_9)
        self.horizontalLayout_18 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_18.setObjectName("horizontalLayout_18")
        self.label_24 = QtWidgets.QLabel(self.horizontalLayoutWidget_19)
        self.label_24.setObjectName("label_24")
        self.horizontalLayout_18.addWidget(self.label_24)
        self.mode_lcms_rb = QtWidgets.QRadioButton(self.horizontalLayoutWidget_19)
        self.mode_lcms_rb.setChecked(True)
        self.mode_lcms_rb.setAutoExclusive(True)
        self.mode_lcms_rb.setObjectName("mode_lcms_rb")
        self.horizontalLayout_18.addWidget(self.mode_lcms_rb)
        self.mode_static_rb = QtWidgets.QRadioButton(self.horizontalLayoutWidget_19)
        self.mode_static_rb.setChecked(False)
        self.mode_static_rb.setObjectName("mode_static_rb")
        self.horizontalLayout_18.addWidget(self.mode_static_rb)
        self.horizontalLayout_14.addLayout(self.horizontalLayout_18)
        self.horizontalLayout_31.addLayout(self.horizontalLayout_14)
        spacerItem17 = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_31.addItem(spacerItem17)
        self.tabframe.addTab(self.hunter_tab, "")
        self.batch_tab = QtWidgets.QWidget()
        self.batch_tab.setObjectName("batch_tab")
        self.gridLayoutWidget_7 = QtWidgets.QWidget(self.batch_tab)
        self.gridLayoutWidget_7.setGeometry(
            QtCore.QRect(10, 10, 901 * scale, 581 * scale)
        )
        self.gridLayoutWidget_7.setObjectName("gridLayoutWidget_7")
        self.gridLayout_9 = QtWidgets.QGridLayout(self.gridLayoutWidget_7)
        self.gridLayout_9.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_9.setObjectName("gridLayout_9")
        self.line_23 = QtWidgets.QFrame(self.gridLayoutWidget_7)
        self.line_23.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_23.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_23.setObjectName("line_23")
        self.gridLayout_9.addWidget(self.line_23, 2, 1, 1, 1)
        self.tab_b_statusrun_pte = QtWidgets.QPlainTextEdit(self.gridLayoutWidget_7)
        self.tab_b_statusrun_pte.setMaximumSize(QtCore.QSize(16777215, 400))
        self.tab_b_statusrun_pte.setObjectName("tab_b_statusrun_pte")
        self.gridLayout_9.addWidget(self.tab_b_statusrun_pte, 6, 1, 1, 1)
        self.label_98 = QtWidgets.QLabel(self.gridLayoutWidget_7)
        self.label_98.setObjectName("label_98")
        self.gridLayout_9.addWidget(self.label_98, 0, 0, 1, 1)
        self.horizontalLayout_26 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_26.setObjectName("horizontalLayout_26")
        self.tab_b_1_lb = QtWidgets.QLabel(self.gridLayoutWidget_7)
        self.tab_b_1_lb.setObjectName("tab_b_1_lb")
        self.horizontalLayout_26.addWidget(self.tab_b_1_lb)
        self.tab_b_addcfg_pb = QtWidgets.QPushButton(self.gridLayoutWidget_7)
        self.tab_b_addcfg_pb.setObjectName("tab_b_addcfg_pb")
        self.horizontalLayout_26.addWidget(self.tab_b_addcfg_pb)
        spacerItem18 = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_26.addItem(spacerItem18)
        self.tab_b_2_lb = QtWidgets.QLabel(self.gridLayoutWidget_7)
        self.tab_b_2_lb.setObjectName("tab_b_2_lb")
        self.horizontalLayout_26.addWidget(self.tab_b_2_lb)
        self.tab_b_addcfgfolder_pb = QtWidgets.QPushButton(self.gridLayoutWidget_7)
        self.tab_b_addcfgfolder_pb.setObjectName("tab_b_addcfgfolder_pb")
        self.horizontalLayout_26.addWidget(self.tab_b_addcfgfolder_pb)
        self.gridLayout_9.addLayout(self.horizontalLayout_26, 0, 1, 1, 1)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        spacerItem19 = QtWidgets.QSpacerItem(
            20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
        )
        self.verticalLayout_2.addItem(spacerItem19)
        self.tab_b_clearall_pb = QtWidgets.QPushButton(self.gridLayoutWidget_7)
        self.tab_b_clearall_pb.setObjectName("tab_b_clearall_pb")
        self.verticalLayout_2.addWidget(self.tab_b_clearall_pb)
        spacerItem20 = QtWidgets.QSpacerItem(
            20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
        )
        self.verticalLayout_2.addItem(spacerItem20)
        self.gridLayout_9.addLayout(self.verticalLayout_2, 1, 2, 1, 1)
        self.line_22 = QtWidgets.QFrame(self.gridLayoutWidget_7)
        self.line_22.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_22.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_22.setObjectName("line_22")
        self.gridLayout_9.addWidget(self.line_22, 5, 1, 1, 1)
        self.verticalLayout_3 = QtWidgets.QVBoxLayout()
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.tab_b_runbatch_pb = QtWidgets.QPushButton(self.gridLayoutWidget_7)
        self.tab_b_runbatch_pb.setObjectName("tab_b_runbatch_pb")
        self.verticalLayout_3.addWidget(self.tab_b_runbatch_pb)
        self.tab_b_runbatch_pgb = QtWidgets.QProgressBar(self.gridLayoutWidget_7)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.tab_b_runbatch_pgb.sizePolicy().hasHeightForWidth()
        )
        self.tab_b_runbatch_pgb.setSizePolicy(sizePolicy)
        self.tab_b_runbatch_pgb.setMaximumSize(QtCore.QSize(160, 16777215))
        self.tab_b_runbatch_pgb.setProperty("value", 5)
        self.tab_b_runbatch_pgb.setAlignment(QtCore.Qt.AlignCenter)
        self.tab_b_runbatch_pgb.setTextVisible(False)
        self.tab_b_runbatch_pgb.setObjectName("tab_b_runbatch_pgb")
        self.verticalLayout_3.addWidget(self.tab_b_runbatch_pgb)
        self.gridLayout_9.addLayout(self.verticalLayout_3, 6, 0, 1, 1)
        self.horizontalLayout_25 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_25.setObjectName("horizontalLayout_25")
        self.tab_b_maxsubcore_lb = QtWidgets.QLabel(self.gridLayoutWidget_7)
        self.tab_b_maxsubcore_lb.setObjectName("tab_b_maxsubcore_lb")
        self.horizontalLayout_25.addWidget(self.tab_b_maxsubcore_lb)
        self.tab_b_maxsubcore_spb = QtWidgets.QSpinBox(self.gridLayoutWidget_7)
        self.tab_b_maxsubcore_spb.setMaximum(16)
        self.tab_b_maxsubcore_spb.setProperty("value", 4)
        self.tab_b_maxsubcore_spb.setObjectName("tab_b_maxsubcore_spb")
        self.horizontalLayout_25.addWidget(self.tab_b_maxsubcore_spb)
        spacerItem21 = QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_25.addItem(spacerItem21)
        self.tab_b_maxsubram_lb = QtWidgets.QLabel(self.gridLayoutWidget_7)
        self.tab_b_maxsubram_lb.setObjectName("tab_b_maxsubram_lb")
        self.horizontalLayout_25.addWidget(self.tab_b_maxsubram_lb)
        self.tab_b_maxsubram_spb = QtWidgets.QSpinBox(self.gridLayoutWidget_7)
        self.tab_b_maxsubram_spb.setMinimum(1)
        self.tab_b_maxsubram_spb.setMaximum(64)
        self.tab_b_maxsubram_spb.setProperty("value", 6)
        self.tab_b_maxsubram_spb.setObjectName("tab_b_maxsubram_spb")
        self.horizontalLayout_25.addWidget(self.tab_b_maxsubram_spb)
        spacerItem22 = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_25.addItem(spacerItem22)
        self.gridLayout_9.addLayout(self.horizontalLayout_25, 3, 1, 1, 1)
        self.horizontalLayout_27 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_27.setObjectName("horizontalLayout_27")
        self.tab_b_mutlimode_cmb = QtWidgets.QComboBox(self.gridLayoutWidget_7)
        self.tab_b_mutlimode_cmb.setMinimumSize(QtCore.QSize(250, 0))
        self.tab_b_mutlimode_cmb.setObjectName("tab_b_mutlimode_cmb")
        self.tab_b_mutlimode_cmb.addItem("")
        self.tab_b_mutlimode_cmb.addItem("")
        self.horizontalLayout_27.addWidget(self.tab_b_mutlimode_cmb)
        self.tab_b_maxbatch_lb = QtWidgets.QLabel(self.gridLayoutWidget_7)
        self.tab_b_maxbatch_lb.setObjectName("tab_b_maxbatch_lb")
        self.horizontalLayout_27.addWidget(self.tab_b_maxbatch_lb)
        self.tab_b_maxbatch_spb = QtWidgets.QSpinBox(self.gridLayoutWidget_7)
        self.tab_b_maxbatch_spb.setMaximum(16)
        self.tab_b_maxbatch_spb.setProperty("value", 3)
        self.tab_b_maxbatch_spb.setObjectName("tab_b_maxbatch_spb")
        self.horizontalLayout_27.addWidget(self.tab_b_maxbatch_spb)
        spacerItem23 = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_27.addItem(spacerItem23)
        self.gridLayout_9.addLayout(self.horizontalLayout_27, 4, 1, 1, 1)
        self.tab_b_infiles_pte = QtWidgets.QPlainTextEdit(self.gridLayoutWidget_7)
        self.tab_b_infiles_pte.setObjectName("tab_b_infiles_pte")
        self.gridLayout_9.addWidget(self.tab_b_infiles_pte, 1, 1, 1, 1)
        self.label_99 = QtWidgets.QLabel(self.gridLayoutWidget_7)
        self.label_99.setObjectName("label_99")
        self.gridLayout_9.addWidget(self.label_99, 3, 0, 1, 1)
        self.tabframe.addTab(self.batch_tab, "")
        self.cfg_tab = QtWidgets.QWidget()
        self.cfg_tab.setObjectName("cfg_tab")
        self.gridLayoutWidget_6 = QtWidgets.QWidget(self.cfg_tab)
        self.gridLayoutWidget_6.setGeometry(
            QtCore.QRect(10, 10, 901 * scale, 320 * scale)
        )
        self.gridLayoutWidget_6.setObjectName("gridLayoutWidget_6")
        self.gridLayout_6 = QtWidgets.QGridLayout(self.gridLayoutWidget_6)
        self.gridLayout_6.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_6.setObjectName("gridLayout_6")
        self.tab_c_hgcfg_le = QtWidgets.QLineEdit(self.gridLayoutWidget_6)
        self.tab_c_hgcfg_le.setObjectName("tab_c_hgcfg_le")
        self.gridLayout_6.addWidget(self.tab_c_hgcfg_le, 5, 1, 1, 1)
        self.horizontalLayout_23 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_23.setObjectName("horizontalLayout_23")
        self.label_93 = QtWidgets.QLabel(self.gridLayoutWidget_6)
        self.label_93.setObjectName("label_93")
        self.horizontalLayout_23.addWidget(self.label_93)
        self.tab_c_imagetype_cmb = QtWidgets.QComboBox(self.gridLayoutWidget_6)
        self.tab_c_imagetype_cmb.setObjectName("tab_c_imagetype_cmb")
        self.tab_c_imagetype_cmb.addItem("")
        self.tab_c_imagetype_cmb.addItem("")
        self.horizontalLayout_23.addWidget(self.tab_c_imagetype_cmb)
        self.label_94 = QtWidgets.QLabel(self.gridLayoutWidget_6)
        self.label_94.setObjectName("label_94")
        self.horizontalLayout_23.addWidget(self.label_94)
        self.tab_c_dpi_spb = QtWidgets.QSpinBox(self.gridLayoutWidget_6)
        self.tab_c_dpi_spb.setMinimum(50)
        self.tab_c_dpi_spb.setMaximum(1000)
        self.tab_c_dpi_spb.setSingleStep(100)
        self.tab_c_dpi_spb.setProperty("value", 300)
        self.tab_c_dpi_spb.setObjectName("tab_c_dpi_spb")
        self.horizontalLayout_23.addWidget(self.tab_c_dpi_spb)
        spacerItem24 = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_23.addItem(spacerItem24)
        self.gridLayout_6.addLayout(self.horizontalLayout_23, 8, 1, 1, 1)
        self.line_12 = QtWidgets.QFrame(self.gridLayoutWidget_6)
        self.line_12.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_12.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_12.setObjectName("line_12")
        self.gridLayout_6.addWidget(self.line_12, 6, 1, 1, 1)
        self.horizontalLayout_36 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_36.setObjectName("horizontalLayout_36")
        self.tab_c_scorecfgtg_lb = QtWidgets.QLabel(self.gridLayoutWidget_6)
        self.tab_c_scorecfgtg_lb.setMinimumSize(QtCore.QSize(80, 0))
        self.tab_c_scorecfgtg_lb.setObjectName("tab_c_scorecfgtg_lb")
        self.horizontalLayout_36.addWidget(self.tab_c_scorecfgtg_lb)
        self.tab_c_scorecfgtg_le = QtWidgets.QLineEdit(self.gridLayoutWidget_6)
        self.tab_c_scorecfgtg_le.setObjectName("tab_c_scorecfgtg_le")
        self.horizontalLayout_36.addWidget(self.tab_c_scorecfgtg_le)
        self.gridLayout_6.addLayout(self.horizontalLayout_36, 3, 1, 1, 1)
        self.horizontalLayout_43 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_43.setObjectName("horizontalLayout_43")
        self.tab_c_scorecfgdg_lb = QtWidgets.QLabel(self.gridLayoutWidget_6)
        self.tab_c_scorecfgdg_lb.setMinimumSize(QtCore.QSize(80, 0))
        self.tab_c_scorecfgdg_lb.setObjectName("tab_c_scorecfgdg_lb")
        self.horizontalLayout_43.addWidget(self.tab_c_scorecfgdg_lb)
        self.tab_c_scorecfgdg_le = QtWidgets.QLineEdit(self.gridLayoutWidget_6)
        self.tab_c_scorecfgdg_le.setObjectName("tab_c_scorecfgdg_le")
        self.horizontalLayout_43.addWidget(self.tab_c_scorecfgdg_le)
        self.gridLayout_6.addLayout(self.horizontalLayout_43, 4, 1, 1, 1)
        self.tab_c_isotopescoremode_lb = QtWidgets.QLabel(self.gridLayoutWidget_6)
        self.tab_c_isotopescoremode_lb.setObjectName("tab_c_isotopescoremode_lb")
        self.gridLayout_6.addWidget(self.tab_c_isotopescoremode_lb, 10, 0, 1, 1)
        self.label_14 = QtWidgets.QLabel(self.gridLayoutWidget_6)
        self.label_14.setObjectName("label_14")
        self.gridLayout_6.addWidget(self.label_14, 5, 0, 1, 1)
        self.tab_c_scorecfgdg_pb = QtWidgets.QPushButton(self.gridLayoutWidget_6)
        self.tab_c_scorecfgdg_pb.setObjectName("tab_c_scorecfgdg_pb")
        self.gridLayout_6.addWidget(self.tab_c_scorecfgdg_pb, 4, 2, 1, 1)
        self.horizontalLayout_37 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_37.setObjectName("horizontalLayout_37")
        self.tab_c_scorecfgpl_lb = QtWidgets.QLabel(self.gridLayoutWidget_6)
        self.tab_c_scorecfgpl_lb.setMinimumSize(QtCore.QSize(80, 0))
        self.tab_c_scorecfgpl_lb.setObjectName("tab_c_scorecfgpl_lb")
        self.horizontalLayout_37.addWidget(self.tab_c_scorecfgpl_lb)
        self.tab_c_scorecfgpl_le = QtWidgets.QLineEdit(self.gridLayoutWidget_6)
        self.tab_c_scorecfgpl_le.setObjectName("tab_c_scorecfgpl_le")
        self.horizontalLayout_37.addWidget(self.tab_c_scorecfgpl_le)
        self.gridLayout_6.addLayout(self.horizontalLayout_37, 2, 1, 1, 1)
        self.tab_c_scorecfgpl_pb = QtWidgets.QPushButton(self.gridLayoutWidget_6)
        self.tab_c_scorecfgpl_pb.setObjectName("tab_c_scorecfgpl_pb")
        self.gridLayout_6.addWidget(self.tab_c_scorecfgpl_pb, 2, 2, 1, 1)
        self.tab_c_scorecfgtg_pb = QtWidgets.QPushButton(self.gridLayoutWidget_6)
        self.tab_c_scorecfgtg_pb.setObjectName("tab_c_scorecfgtg_pb")
        self.gridLayout_6.addWidget(self.tab_c_scorecfgtg_pb, 3, 2, 1, 1)
        self.tab_c_hgcfg_pb = QtWidgets.QPushButton(self.gridLayoutWidget_6)
        self.tab_c_hgcfg_pb.setObjectName("tab_c_hgcfg_pb")
        self.gridLayout_6.addWidget(self.tab_c_hgcfg_pb, 5, 2, 1, 1)
        self.label_21 = QtWidgets.QLabel(self.gridLayoutWidget_6)
        self.label_21.setObjectName("label_21")
        self.gridLayout_6.addWidget(self.label_21, 2, 0, 1, 1)
        self.horizontalLayout_22 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_22.setObjectName("horizontalLayout_22")
        self.tab_c_parallization_cmb = QtWidgets.QComboBox(self.gridLayoutWidget_6)
        self.tab_c_parallization_cmb.setMinimumSize(QtCore.QSize(200, 0))
        self.tab_c_parallization_cmb.setObjectName("tab_c_parallization_cmb")
        self.tab_c_parallization_cmb.addItem("")
        self.tab_c_parallization_cmb.addItem("")
        self.horizontalLayout_22.addWidget(self.tab_c_parallization_cmb)
        self.label_63 = QtWidgets.QLabel(self.gridLayoutWidget_6)
        self.label_63.setObjectName("label_63")
        self.horizontalLayout_22.addWidget(self.label_63)
        self.tab_c_cores_spb = QtWidgets.QSpinBox(self.gridLayoutWidget_6)
        self.tab_c_cores_spb.setMinimum(1)
        self.tab_c_cores_spb.setProperty("value", 3)
        self.tab_c_cores_spb.setObjectName("tab_c_cores_spb")
        self.horizontalLayout_22.addWidget(self.tab_c_cores_spb)
        self.label_64 = QtWidgets.QLabel(self.gridLayoutWidget_6)
        self.label_64.setObjectName("label_64")
        self.horizontalLayout_22.addWidget(self.label_64)
        self.tab_c_ram_spb = QtWidgets.QSpinBox(self.gridLayoutWidget_6)
        self.tab_c_ram_spb.setMinimum(5)
        self.tab_c_ram_spb.setProperty("value", 5)
        self.tab_c_ram_spb.setObjectName("tab_c_ram_spb")
        self.horizontalLayout_22.addWidget(self.tab_c_ram_spb)
        spacerItem25 = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_22.addItem(spacerItem25)
        self.gridLayout_6.addLayout(self.horizontalLayout_22, 7, 1, 1, 1)
        self.horizontalLayout_30 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_30.setObjectName("horizontalLayout_30")
        self.tab_c_tag_all_fa_chb = QtWidgets.QCheckBox(self.gridLayoutWidget_6)
        self.tab_c_tag_all_fa_chb.setChecked(True)
        self.tab_c_tag_all_fa_chb.setObjectName("tab_c_tag_all_fa_chb")
        self.horizontalLayout_30.addWidget(self.tab_c_tag_all_fa_chb)
        spacerItem26 = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_30.addItem(spacerItem26)
        self.gridLayout_6.addLayout(self.horizontalLayout_30, 9, 1, 1, 1)
        self.label_12 = QtWidgets.QLabel(self.gridLayoutWidget_6)
        self.label_12.setObjectName("label_12")
        self.gridLayout_6.addWidget(self.label_12, 0, 0, 1, 1)
        self.tab_c_falistpl_pb = QtWidgets.QPushButton(self.gridLayoutWidget_6)
        self.tab_c_falistpl_pb.setObjectName("tab_c_falistpl_pb")
        self.gridLayout_6.addWidget(self.tab_c_falistpl_pb, 0, 2, 1, 1)
        self.tab_c_savesettings_pb = QtWidgets.QPushButton(self.gridLayoutWidget_6)
        self.tab_c_savesettings_pb.setObjectName("tab_c_savesettings_pb")
        self.gridLayout_6.addWidget(self.tab_c_savesettings_pb, 11, 1, 1, 1)
        spacerItem27 = QtWidgets.QSpacerItem(
            190, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum
        )
        self.gridLayout_6.addItem(spacerItem27, 11, 0, 1, 1)
        self.horizontalLayout_38 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_38.setObjectName("horizontalLayout_38")
        self.tab_c_falisttg_lb = QtWidgets.QLabel(self.gridLayoutWidget_6)
        self.tab_c_falisttg_lb.setMinimumSize(QtCore.QSize(80, 0))
        self.tab_c_falisttg_lb.setObjectName("tab_c_falisttg_lb")
        self.horizontalLayout_38.addWidget(self.tab_c_falisttg_lb)
        self.tab_c_falisttg_le = QtWidgets.QLineEdit(self.gridLayoutWidget_6)
        self.tab_c_falisttg_le.setObjectName("tab_c_falisttg_le")
        self.horizontalLayout_38.addWidget(self.tab_c_falisttg_le)
        self.gridLayout_6.addLayout(self.horizontalLayout_38, 1, 1, 1, 1)
        self.horizontalLayout_39 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_39.setObjectName("horizontalLayout_39")
        self.tab_c_falistpl_lb = QtWidgets.QLabel(self.gridLayoutWidget_6)
        self.tab_c_falistpl_lb.setMinimumSize(QtCore.QSize(80, 0))
        self.tab_c_falistpl_lb.setObjectName("tab_c_falistpl_lb")
        self.horizontalLayout_39.addWidget(self.tab_c_falistpl_lb)
        self.tab_c_falistpl_le = QtWidgets.QLineEdit(self.gridLayoutWidget_6)
        self.tab_c_falistpl_le.setObjectName("tab_c_falistpl_le")
        self.horizontalLayout_39.addWidget(self.tab_c_falistpl_le)
        self.gridLayout_6.addLayout(self.horizontalLayout_39, 0, 1, 1, 1)
        self.tab_c_falisttg_pb = QtWidgets.QPushButton(self.gridLayoutWidget_6)
        self.tab_c_falisttg_pb.setObjectName("tab_c_falisttg_pb")
        self.gridLayout_6.addWidget(self.tab_c_falisttg_pb, 1, 2, 1, 1)
        self.horizontalLayout_24 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_24.setObjectName("horizontalLayout_24")
        self.tab_c_isotopescoremode_cmb = QtWidgets.QComboBox(self.gridLayoutWidget_6)
        self.tab_c_isotopescoremode_cmb.setMinimumSize(QtCore.QSize(150, 0))
        self.tab_c_isotopescoremode_cmb.setObjectName("tab_c_isotopescoremode_cmb")
        self.tab_c_isotopescoremode_cmb.addItem("")
        self.tab_c_isotopescoremode_cmb.addItem("")
        self.horizontalLayout_24.addWidget(self.tab_c_isotopescoremode_cmb)
        self.tab_c_scoremode_lb = QtWidgets.QLabel(self.gridLayoutWidget_6)
        self.tab_c_scoremode_lb.setObjectName("tab_c_scoremode_lb")
        self.horizontalLayout_24.addWidget(self.tab_c_scoremode_lb)
        self.tab_c_scoremode_cmb = QtWidgets.QComboBox(self.gridLayoutWidget_6)
        self.tab_c_scoremode_cmb.setMinimumSize(QtCore.QSize(150, 0))
        self.tab_c_scoremode_cmb.setObjectName("tab_c_scoremode_cmb")
        self.tab_c_scoremode_cmb.addItem("")
        self.tab_c_scoremode_cmb.addItem("")
        self.horizontalLayout_24.addWidget(self.tab_c_scoremode_cmb)
        spacerItem28 = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_24.addItem(spacerItem28)
        self.gridLayout_6.addLayout(self.horizontalLayout_24, 10, 1, 1, 1)
        self.label_35 = QtWidgets.QLabel(self.gridLayoutWidget_6)
        self.label_35.setObjectName("label_35")
        self.gridLayout_6.addWidget(self.label_35, 7, 0, 1, 1)
        self.label_54 = QtWidgets.QLabel(self.gridLayoutWidget_6)
        self.label_54.setObjectName("label_54")
        self.gridLayout_6.addWidget(self.label_54, 8, 0, 1, 1)
        self.gridLayoutWidget = QtWidgets.QWidget(self.cfg_tab)
        self.gridLayoutWidget.setGeometry(
            QtCore.QRect(10, 360 * scale, 901 * scale, 251 * scale)
        )
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout_2.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.tab_c_lmcalcfalist_le = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.tab_c_lmcalcfalist_le.setObjectName("tab_c_lmcalcfalist_le")
        self.gridLayout_2.addWidget(self.tab_c_lmcalcfalist_le, 2, 1, 1, 1)
        self.tab_c_lmexport_le = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.tab_c_lmexport_le.setObjectName("tab_c_lmexport_le")
        self.gridLayout_2.addWidget(self.tab_c_lmexport_le, 3, 1, 1, 1)
        self.tab_c_lmstatus_pte = QtWidgets.QPlainTextEdit(self.gridLayoutWidget)
        self.tab_c_lmstatus_pte.setObjectName("tab_c_lmstatus_pte")
        self.gridLayout_2.addWidget(self.tab_c_lmstatus_pte, 4, 1, 1, 1)
        self.label_28 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_28.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.label_28.setObjectName("label_28")
        self.gridLayout_2.addWidget(self.label_28, 2, 0, 1, 1)
        self.tab_c_lmcalcfalist_pb = QtWidgets.QPushButton(self.gridLayoutWidget)
        self.tab_c_lmcalcfalist_pb.setObjectName("tab_c_lmcalcfalist_pb")
        self.gridLayout_2.addWidget(self.tab_c_lmcalcfalist_pb, 2, 2, 1, 1)
        self.tab_c_lmexport_pb = QtWidgets.QPushButton(self.gridLayoutWidget)
        self.tab_c_lmexport_pb.setObjectName("tab_c_lmexport_pb")
        self.gridLayout_2.addWidget(self.tab_c_lmexport_pb, 3, 2, 1, 1)
        self.label_29 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_29.setObjectName("label_29")
        self.gridLayout_2.addWidget(self.label_29, 0, 1, 1, 1)
        self.label_27 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_27.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.label_27.setObjectName("label_27")
        self.gridLayout_2.addWidget(self.label_27, 3, 0, 1, 1)
        self.label_48 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_48.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.label_48.setObjectName("label_48")
        self.gridLayout_2.addWidget(self.label_48, 1, 0, 1, 1)
        self.horizontalLayout_44 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_44.setObjectName("horizontalLayout_44")
        self.tab_c_lipidclass_cmb = QtWidgets.QComboBox(self.gridLayoutWidget)
        self.tab_c_lipidclass_cmb.setMinimumSize(QtCore.QSize(220, 0))
        self.tab_c_lipidclass_cmb.setMaxVisibleItems(15)
        self.tab_c_lipidclass_cmb.setObjectName("tab_c_lipidclass_cmb")
        self.tab_c_lipidclass_cmb.addItem("")
        self.tab_c_lipidclass_cmb.addItem("")
        self.tab_c_lipidclass_cmb.addItem("")
        self.tab_c_lipidclass_cmb.addItem("")
        self.tab_c_lipidclass_cmb.addItem("")
        self.tab_c_lipidclass_cmb.addItem("")
        self.tab_c_lipidclass_cmb.addItem("")
        self.tab_c_lipidclass_cmb.addItem("")
        self.tab_c_lipidclass_cmb.addItem("")
        self.tab_c_lipidclass_cmb.addItem("")
        self.tab_c_lipidclass_cmb.addItem("")
        self.tab_c_lipidclass_cmb.addItem("")
        self.tab_c_lipidclass_cmb.addItem("")
        self.tab_c_lipidclass_cmb.addItem("")
        self.tab_c_lipidclass_cmb.addItem("")
        self.tab_c_lipidclass_cmb.addItem("")
        self.tab_c_lipidclass_cmb.addItem("")
        self.tab_c_lipidclass_cmb.addItem("")
        self.tab_c_lipidclass_cmb.addItem("")
        self.horizontalLayout_44.addWidget(self.tab_c_lipidclass_cmb)
        self.label_17 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_17.setObjectName("label_17")
        self.horizontalLayout_44.addWidget(self.label_17)
        self.tab_c_lmms2ppm_spb = QtWidgets.QSpinBox(self.gridLayoutWidget)
        self.tab_c_lmms2ppm_spb.setMaximum(999)
        self.tab_c_lmms2ppm_spb.setProperty("value", 50)
        self.tab_c_lmms2ppm_spb.setObjectName("tab_c_lmms2ppm_spb")
        self.horizontalLayout_44.addWidget(self.tab_c_lmms2ppm_spb)
        spacerItem29 = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_44.addItem(spacerItem29)
        self.gridLayout_2.addLayout(self.horizontalLayout_44, 1, 1, 1, 1)
        self.verticalLayout_6 = QtWidgets.QVBoxLayout()
        self.verticalLayout_6.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.tab_c_lmrun_pb = QtWidgets.QPushButton(self.gridLayoutWidget)
        self.tab_c_lmrun_pb.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.tab_c_lmrun_pb.setObjectName("tab_c_lmrun_pb")
        self.verticalLayout_6.addWidget(self.tab_c_lmrun_pb)
        self.tab_c_runlm_pgb = QtWidgets.QProgressBar(self.gridLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.tab_c_runlm_pgb.sizePolicy().hasHeightForWidth()
        )
        self.tab_c_runlm_pgb.setSizePolicy(sizePolicy)
        self.tab_c_runlm_pgb.setMaximumSize(QtCore.QSize(190, 16777215))
        self.tab_c_runlm_pgb.setProperty("value", 5)
        self.tab_c_runlm_pgb.setAlignment(QtCore.Qt.AlignCenter)
        self.tab_c_runlm_pgb.setTextVisible(False)
        self.tab_c_runlm_pgb.setObjectName("tab_c_runlm_pgb")
        self.verticalLayout_6.addWidget(self.tab_c_runlm_pgb)
        self.gridLayout_2.addLayout(self.verticalLayout_6, 4, 0, 1, 1)
        self.line = QtWidgets.QFrame(self.cfg_tab)
        self.line.setGeometry(QtCore.QRect(10, 340 * scale, 901 * scale, 20 * scale))
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.tabframe.addTab(self.cfg_tab, "")
        self.generator_tab = QtWidgets.QWidget()
        self.generator_tab.setObjectName("generator_tab")
        self.lipidgen_tabframe = QtWidgets.QTabWidget(self.generator_tab)
        self.lipidgen_tabframe.setGeometry(
            QtCore.QRect(10, 0, 911 * scale, 611 * scale)
        )
        self.lipidgen_tabframe.setObjectName("lipidgen_tabframe")
        self.lipidgen_tab1 = QtWidgets.QWidget()
        self.lipidgen_tab1.setObjectName("lipidgen_tab1")
        self.horizontalLayoutWidget_7 = QtWidgets.QWidget(self.lipidgen_tab1)
        self.horizontalLayoutWidget_7.setGeometry(
            QtCore.QRect(20, 0, 861 * scale, 31 * scale)
        )
        self.horizontalLayoutWidget_7.setObjectName("horizontalLayoutWidget_7")
        self.horizontalLayout_5 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_7)
        self.horizontalLayout_5.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.label_2 = QtWidgets.QLabel(self.horizontalLayoutWidget_7)
        self.label_2.setObjectName("label_2")
        self.horizontalLayout_5.addWidget(self.label_2)
        self.lipidgen_tab1_lipidclass_cmb = QtWidgets.QComboBox(
            self.horizontalLayoutWidget_7
        )
        self.lipidgen_tab1_lipidclass_cmb.setMinimumSize(QtCore.QSize(250, 0))
        self.lipidgen_tab1_lipidclass_cmb.setObjectName("lipidgen_tab1_lipidclass_cmb")
        self.lipidgen_tab1_lipidclass_cmb.addItem("")
        self.lipidgen_tab1_lipidclass_cmb.addItem("")
        self.lipidgen_tab1_lipidclass_cmb.addItem("")
        self.lipidgen_tab1_lipidclass_cmb.addItem("")
        self.lipidgen_tab1_lipidclass_cmb.addItem("")
        self.lipidgen_tab1_lipidclass_cmb.addItem("")
        self.horizontalLayout_5.addWidget(self.lipidgen_tab1_lipidclass_cmb)
        self.lipidgen_tab1_all_sn_chb = QtWidgets.QCheckBox(
            self.horizontalLayoutWidget_7
        )
        self.lipidgen_tab1_all_sn_chb.setObjectName("lipidgen_tab1_all_sn_chb")
        self.horizontalLayout_5.addWidget(self.lipidgen_tab1_all_sn_chb)
        spacerItem30 = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_5.addItem(spacerItem30)
        self.horizontalLayoutWidget_9 = QtWidgets.QWidget(self.lipidgen_tab1)
        self.horizontalLayoutWidget_9.setGeometry(
            QtCore.QRect(20, 510 * scale, 861 * scale, 31 * scale)
        )
        self.horizontalLayoutWidget_9.setObjectName("horizontalLayoutWidget_9")
        self.horizontalLayout_8 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_9)
        self.horizontalLayout_8.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_8.setObjectName("horizontalLayout_8")
        self.label_59 = QtWidgets.QLabel(self.horizontalLayoutWidget_9)
        self.label_59.setMinimumSize(QtCore.QSize(120, 0))
        self.label_59.setObjectName("label_59")
        self.horizontalLayout_8.addWidget(self.label_59)
        self.lipidgen_tab1_savelist_le = QtWidgets.QLineEdit(
            self.horizontalLayoutWidget_9
        )
        self.lipidgen_tab1_savelist_le.setObjectName("lipidgen_tab1_savelist_le")
        self.horizontalLayout_8.addWidget(self.lipidgen_tab1_savelist_le)
        self.lipidgen_tab1_savelist_pb = QtWidgets.QPushButton(
            self.horizontalLayoutWidget_9
        )
        self.lipidgen_tab1_savelist_pb.setObjectName("lipidgen_tab1_savelist_pb")
        self.horizontalLayout_8.addWidget(self.lipidgen_tab1_savelist_pb)
        self.verticalLayoutWidget_3 = QtWidgets.QWidget(self.lipidgen_tab1)
        self.verticalLayoutWidget_3.setGeometry(
            QtCore.QRect(20, 30 * scale, 861 * scale, 146 * scale)
        )
        self.verticalLayoutWidget_3.setObjectName("verticalLayoutWidget_3")
        self.lipidgen_tab1_sn1_verticalLayout = QtWidgets.QVBoxLayout(
            self.verticalLayoutWidget_3
        )
        self.lipidgen_tab1_sn1_verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.lipidgen_tab1_sn1_verticalLayout.setObjectName(
            "lipidgen_tab1_sn1_verticalLayout"
        )
        self.horizontalLayout_19 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_19.setObjectName("horizontalLayout_19")
        self.lipidgen_tab1_sn1_lb = QtWidgets.QLabel(self.verticalLayoutWidget_3)
        self.lipidgen_tab1_sn1_lb.setObjectName("lipidgen_tab1_sn1_lb")
        self.horizontalLayout_19.addWidget(self.lipidgen_tab1_sn1_lb)
        spacerItem31 = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_19.addItem(spacerItem31)
        self.lipidgen_tab1_sn1allfa_chb = QtWidgets.QCheckBox(
            self.verticalLayoutWidget_3
        )
        self.lipidgen_tab1_sn1allfa_chb.setObjectName("lipidgen_tab1_sn1allfa_chb")
        self.horizontalLayout_19.addWidget(self.lipidgen_tab1_sn1allfa_chb)
        self.lipidgen_tab1_sn1_verticalLayout.addLayout(self.horizontalLayout_19)
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.sn1_fa12x0_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_3)
        self.sn1_fa12x0_chb.setObjectName("sn1_fa12x0_chb")
        self.gridLayout.addWidget(self.sn1_fa12x0_chb, 1, 0, 1, 1)
        self.sn1_fa22x5_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_3)
        self.sn1_fa22x5_chb.setObjectName("sn1_fa22x5_chb")
        self.gridLayout.addWidget(self.sn1_fa22x5_chb, 2, 5, 1, 1)
        self.sn1_fa18x2_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_3)
        self.sn1_fa18x2_chb.setChecked(True)
        self.sn1_fa18x2_chb.setObjectName("sn1_fa18x2_chb")
        self.gridLayout.addWidget(self.sn1_fa18x2_chb, 1, 6, 1, 1)
        self.sn1_fa18x1_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_3)
        self.sn1_fa18x1_chb.setChecked(True)
        self.sn1_fa18x1_chb.setObjectName("sn1_fa18x1_chb")
        self.gridLayout.addWidget(self.sn1_fa18x1_chb, 1, 5, 1, 1)
        self.sn1_fa22x4_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_3)
        self.sn1_fa22x4_chb.setObjectName("sn1_fa22x4_chb")
        self.gridLayout.addWidget(self.sn1_fa22x4_chb, 2, 4, 1, 1)
        self.sn1_fa20x5_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_3)
        self.sn1_fa20x5_chb.setObjectName("sn1_fa20x5_chb")
        self.gridLayout.addWidget(self.sn1_fa20x5_chb, 2, 3, 1, 1)
        self.sn1_fa18x0_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_3)
        self.sn1_fa18x0_chb.setChecked(True)
        self.sn1_fa18x0_chb.setObjectName("sn1_fa18x0_chb")
        self.gridLayout.addWidget(self.sn1_fa18x0_chb, 1, 4, 1, 1)
        self.sn1_fa22x6_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_3)
        self.sn1_fa22x6_chb.setObjectName("sn1_fa22x6_chb")
        self.gridLayout.addWidget(self.sn1_fa22x6_chb, 2, 6, 1, 1)
        self.sn1_fa16x1_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_3)
        self.sn1_fa16x1_chb.setObjectName("sn1_fa16x1_chb")
        self.gridLayout.addWidget(self.sn1_fa16x1_chb, 1, 3, 1, 1)
        self.sn1_fa14x0_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_3)
        self.sn1_fa14x0_chb.setObjectName("sn1_fa14x0_chb")
        self.gridLayout.addWidget(self.sn1_fa14x0_chb, 1, 1, 1, 1)
        self.sn1_fa20x3_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_3)
        self.sn1_fa20x3_chb.setObjectName("sn1_fa20x3_chb")
        self.gridLayout.addWidget(self.sn1_fa20x3_chb, 2, 1, 1, 1)
        self.sn1_fa16x0_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_3)
        self.sn1_fa16x0_chb.setChecked(True)
        self.sn1_fa16x0_chb.setObjectName("sn1_fa16x0_chb")
        self.gridLayout.addWidget(self.sn1_fa16x0_chb, 1, 2, 1, 1)
        self.sn1_fa18x3_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_3)
        self.sn1_fa18x3_chb.setObjectName("sn1_fa18x3_chb")
        self.gridLayout.addWidget(self.sn1_fa18x3_chb, 2, 0, 1, 1)
        self.sn1_fa20x4_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_3)
        self.sn1_fa20x4_chb.setObjectName("sn1_fa20x4_chb")
        self.gridLayout.addWidget(self.sn1_fa20x4_chb, 2, 2, 1, 1)
        self.lipidgen_tab1_sn1_verticalLayout.addLayout(self.gridLayout)
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.sn1_showop_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_3)
        self.sn1_showop_chb.setMaximumSize(QtCore.QSize(180, 16777215))
        self.sn1_showop_chb.setObjectName("sn1_showop_chb")
        self.horizontalLayout_4.addWidget(self.sn1_showop_chb)
        self.horizontalLayout_6 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.sn1_o16x0_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_3)
        self.sn1_o16x0_chb.setObjectName("sn1_o16x0_chb")
        self.horizontalLayout_6.addWidget(self.sn1_o16x0_chb)
        self.sn1_o18x0_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_3)
        self.sn1_o18x0_chb.setObjectName("sn1_o18x0_chb")
        self.horizontalLayout_6.addWidget(self.sn1_o18x0_chb)
        self.sn1_o20x0_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_3)
        self.sn1_o20x0_chb.setObjectName("sn1_o20x0_chb")
        self.horizontalLayout_6.addWidget(self.sn1_o20x0_chb)
        self.sn1_p16x0_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_3)
        self.sn1_p16x0_chb.setObjectName("sn1_p16x0_chb")
        self.horizontalLayout_6.addWidget(self.sn1_p16x0_chb)
        self.sn1_p18x0_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_3)
        self.sn1_p18x0_chb.setObjectName("sn1_p18x0_chb")
        self.horizontalLayout_6.addWidget(self.sn1_p18x0_chb)
        self.sn1_p20x0_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_3)
        self.sn1_p20x0_chb.setObjectName("sn1_p20x0_chb")
        self.horizontalLayout_6.addWidget(self.sn1_p20x0_chb)
        self.horizontalLayout_4.addLayout(self.horizontalLayout_6)
        self.lipidgen_tab1_sn1_verticalLayout.addLayout(self.horizontalLayout_4)
        self.horizontalLayout_10 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_10.setObjectName("horizontalLayout_10")
        self.label_22 = QtWidgets.QLabel(self.verticalLayoutWidget_3)
        self.label_22.setObjectName("label_22")
        self.horizontalLayout_10.addWidget(self.label_22)
        self.sn1_otherfa_lb = QtWidgets.QLineEdit(self.verticalLayoutWidget_3)
        self.sn1_otherfa_lb.setObjectName("sn1_otherfa_lb")
        self.horizontalLayout_10.addWidget(self.sn1_otherfa_lb)
        self.lipidgen_tab1_sn1_verticalLayout.addLayout(self.horizontalLayout_10)
        self.label_62 = QtWidgets.QLabel(self.verticalLayoutWidget_3)
        self.label_62.setObjectName("label_62")
        self.lipidgen_tab1_sn1_verticalLayout.addWidget(self.label_62)
        self.verticalLayoutWidget_4 = QtWidgets.QWidget(self.lipidgen_tab1)
        self.verticalLayoutWidget_4.setGeometry(
            QtCore.QRect(20, 180 * scale, 861 * scale, 100 * scale)
        )
        self.verticalLayoutWidget_4.setObjectName("verticalLayoutWidget_4")
        self.lipidgen_tab1_sn2_verticalLayout = QtWidgets.QVBoxLayout(
            self.verticalLayoutWidget_4
        )
        self.lipidgen_tab1_sn2_verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.lipidgen_tab1_sn2_verticalLayout.setObjectName(
            "lipidgen_tab1_sn2_verticalLayout"
        )
        self.horizontalLayout_20 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_20.setObjectName("horizontalLayout_20")
        self.lipidgen_tab1_sn2_lb = QtWidgets.QLabel(self.verticalLayoutWidget_4)
        self.lipidgen_tab1_sn2_lb.setObjectName("lipidgen_tab1_sn2_lb")
        self.horizontalLayout_20.addWidget(self.lipidgen_tab1_sn2_lb)
        spacerItem32 = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_20.addItem(spacerItem32)
        self.lipidgen_tab1_sn2allfa_chb = QtWidgets.QCheckBox(
            self.verticalLayoutWidget_4
        )
        self.lipidgen_tab1_sn2allfa_chb.setObjectName("lipidgen_tab1_sn2allfa_chb")
        self.horizontalLayout_20.addWidget(self.lipidgen_tab1_sn2allfa_chb)
        self.lipidgen_tab1_sn2_verticalLayout.addLayout(self.horizontalLayout_20)
        self.gridLayout_7 = QtWidgets.QGridLayout()
        self.gridLayout_7.setObjectName("gridLayout_7")
        self.sn2_fa22x5_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_4)
        self.sn2_fa22x5_chb.setChecked(True)
        self.sn2_fa22x5_chb.setObjectName("sn2_fa22x5_chb")
        self.gridLayout_7.addWidget(self.sn2_fa22x5_chb, 1, 5, 1, 1)
        self.sn2_fa22x4_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_4)
        self.sn2_fa22x4_chb.setChecked(True)
        self.sn2_fa22x4_chb.setObjectName("sn2_fa22x4_chb")
        self.gridLayout_7.addWidget(self.sn2_fa22x4_chb, 1, 4, 1, 1)
        self.sn2_fa20x5_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_4)
        self.sn2_fa20x5_chb.setChecked(True)
        self.sn2_fa20x5_chb.setObjectName("sn2_fa20x5_chb")
        self.gridLayout_7.addWidget(self.sn2_fa20x5_chb, 1, 3, 1, 1)
        self.sn2_fa22x6_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_4)
        self.sn2_fa22x6_chb.setChecked(True)
        self.sn2_fa22x6_chb.setObjectName("sn2_fa22x6_chb")
        self.gridLayout_7.addWidget(self.sn2_fa22x6_chb, 1, 6, 1, 1)
        self.sn2_fa18x1_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_4)
        self.sn2_fa18x1_chb.setChecked(True)
        self.sn2_fa18x1_chb.setObjectName("sn2_fa18x1_chb")
        self.gridLayout_7.addWidget(self.sn2_fa18x1_chb, 0, 5, 1, 1)
        self.sn2_fa14x0_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_4)
        self.sn2_fa14x0_chb.setObjectName("sn2_fa14x0_chb")
        self.gridLayout_7.addWidget(self.sn2_fa14x0_chb, 0, 1, 1, 1)
        self.sn2_fa16x0_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_4)
        self.sn2_fa16x0_chb.setChecked(True)
        self.sn2_fa16x0_chb.setObjectName("sn2_fa16x0_chb")
        self.gridLayout_7.addWidget(self.sn2_fa16x0_chb, 0, 2, 1, 1)
        self.sn2_fa18x3_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_4)
        self.sn2_fa18x3_chb.setChecked(True)
        self.sn2_fa18x3_chb.setObjectName("sn2_fa18x3_chb")
        self.gridLayout_7.addWidget(self.sn2_fa18x3_chb, 1, 0, 1, 1)
        self.sn2_fa18x2_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_4)
        self.sn2_fa18x2_chb.setChecked(True)
        self.sn2_fa18x2_chb.setObjectName("sn2_fa18x2_chb")
        self.gridLayout_7.addWidget(self.sn2_fa18x2_chb, 0, 6, 1, 1)
        self.sn2_fa18x0_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_4)
        self.sn2_fa18x0_chb.setChecked(True)
        self.sn2_fa18x0_chb.setObjectName("sn2_fa18x0_chb")
        self.gridLayout_7.addWidget(self.sn2_fa18x0_chb, 0, 4, 1, 1)
        self.sn2_fa16x1_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_4)
        self.sn2_fa16x1_chb.setObjectName("sn2_fa16x1_chb")
        self.gridLayout_7.addWidget(self.sn2_fa16x1_chb, 0, 3, 1, 1)
        self.sn2_fa20x4_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_4)
        self.sn2_fa20x4_chb.setChecked(True)
        self.sn2_fa20x4_chb.setObjectName("sn2_fa20x4_chb")
        self.gridLayout_7.addWidget(self.sn2_fa20x4_chb, 1, 2, 1, 1)
        self.sn2_fa20x3_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_4)
        self.sn2_fa20x3_chb.setChecked(True)
        self.sn2_fa20x3_chb.setObjectName("sn2_fa20x3_chb")
        self.gridLayout_7.addWidget(self.sn2_fa20x3_chb, 1, 1, 1, 1)
        self.sn2_fa12x0_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_4)
        self.sn2_fa12x0_chb.setObjectName("sn2_fa12x0_chb")
        self.gridLayout_7.addWidget(self.sn2_fa12x0_chb, 0, 0, 1, 1)
        self.lipidgen_tab1_sn2_verticalLayout.addLayout(self.gridLayout_7)
        self.horizontalLayout_32 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_32.setObjectName("horizontalLayout_32")
        self.label_41 = QtWidgets.QLabel(self.verticalLayoutWidget_4)
        self.label_41.setObjectName("label_41")
        self.horizontalLayout_32.addWidget(self.label_41)
        self.sn2_otherfa_lb = QtWidgets.QLineEdit(self.verticalLayoutWidget_4)
        self.sn2_otherfa_lb.setObjectName("sn2_otherfa_lb")
        self.horizontalLayout_32.addWidget(self.sn2_otherfa_lb)
        self.lipidgen_tab1_sn2_verticalLayout.addLayout(self.horizontalLayout_32)
        self.verticalLayoutWidget_5 = QtWidgets.QWidget(self.lipidgen_tab1)
        self.verticalLayoutWidget_5.setGeometry(
            QtCore.QRect(20, 290 * scale, 861 * scale, 101 * scale)
        )
        self.verticalLayoutWidget_5.setObjectName("verticalLayoutWidget_5")
        self.lipidgen_tab1_sn3_verticalLayout = QtWidgets.QVBoxLayout(
            self.verticalLayoutWidget_5
        )
        self.lipidgen_tab1_sn3_verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.lipidgen_tab1_sn3_verticalLayout.setObjectName(
            "lipidgen_tab1_sn3_verticalLayout"
        )
        self.horizontalLayout_21 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_21.setObjectName("horizontalLayout_21")
        self.lipidgen_tab1_sn3_lb = QtWidgets.QLabel(self.verticalLayoutWidget_5)
        self.lipidgen_tab1_sn3_lb.setObjectName("lipidgen_tab1_sn3_lb")
        self.horizontalLayout_21.addWidget(self.lipidgen_tab1_sn3_lb)
        spacerItem33 = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_21.addItem(spacerItem33)
        self.lipidgen_tab1_sn3allfa_chb = QtWidgets.QCheckBox(
            self.verticalLayoutWidget_5
        )
        self.lipidgen_tab1_sn3allfa_chb.setObjectName("lipidgen_tab1_sn3allfa_chb")
        self.horizontalLayout_21.addWidget(self.lipidgen_tab1_sn3allfa_chb)
        self.lipidgen_tab1_sn3_verticalLayout.addLayout(self.horizontalLayout_21)
        self.gridLayout_8 = QtWidgets.QGridLayout()
        self.gridLayout_8.setObjectName("gridLayout_8")
        self.sn3_fa14x0_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_5)
        self.sn3_fa14x0_chb.setObjectName("sn3_fa14x0_chb")
        self.gridLayout_8.addWidget(self.sn3_fa14x0_chb, 0, 1, 1, 1)
        self.sn3_fa22x4_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_5)
        self.sn3_fa22x4_chb.setChecked(True)
        self.sn3_fa22x4_chb.setObjectName("sn3_fa22x4_chb")
        self.gridLayout_8.addWidget(self.sn3_fa22x4_chb, 1, 4, 1, 1)
        self.sn3_fa18x1_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_5)
        self.sn3_fa18x1_chb.setChecked(True)
        self.sn3_fa18x1_chb.setObjectName("sn3_fa18x1_chb")
        self.gridLayout_8.addWidget(self.sn3_fa18x1_chb, 0, 5, 1, 1)
        self.sn3_fa22x5_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_5)
        self.sn3_fa22x5_chb.setChecked(True)
        self.sn3_fa22x5_chb.setObjectName("sn3_fa22x5_chb")
        self.gridLayout_8.addWidget(self.sn3_fa22x5_chb, 1, 5, 1, 1)
        self.sn3_fa22x6_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_5)
        self.sn3_fa22x6_chb.setChecked(True)
        self.sn3_fa22x6_chb.setObjectName("sn3_fa22x6_chb")
        self.gridLayout_8.addWidget(self.sn3_fa22x6_chb, 1, 6, 1, 1)
        self.sn3_fa18x3_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_5)
        self.sn3_fa18x3_chb.setChecked(True)
        self.sn3_fa18x3_chb.setObjectName("sn3_fa18x3_chb")
        self.gridLayout_8.addWidget(self.sn3_fa18x3_chb, 1, 0, 1, 1)
        self.sn3_fa16x1_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_5)
        self.sn3_fa16x1_chb.setChecked(False)
        self.sn3_fa16x1_chb.setObjectName("sn3_fa16x1_chb")
        self.gridLayout_8.addWidget(self.sn3_fa16x1_chb, 0, 3, 1, 1)
        self.sn3_fa20x4_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_5)
        self.sn3_fa20x4_chb.setChecked(True)
        self.sn3_fa20x4_chb.setObjectName("sn3_fa20x4_chb")
        self.gridLayout_8.addWidget(self.sn3_fa20x4_chb, 1, 2, 1, 1)
        self.sn3_fa18x0_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_5)
        self.sn3_fa18x0_chb.setChecked(True)
        self.sn3_fa18x0_chb.setObjectName("sn3_fa18x0_chb")
        self.gridLayout_8.addWidget(self.sn3_fa18x0_chb, 0, 4, 1, 1)
        self.sn3_fa20x5_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_5)
        self.sn3_fa20x5_chb.setChecked(True)
        self.sn3_fa20x5_chb.setObjectName("sn3_fa20x5_chb")
        self.gridLayout_8.addWidget(self.sn3_fa20x5_chb, 1, 3, 1, 1)
        self.sn3_fa18x2_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_5)
        self.sn3_fa18x2_chb.setChecked(True)
        self.sn3_fa18x2_chb.setObjectName("sn3_fa18x2_chb")
        self.gridLayout_8.addWidget(self.sn3_fa18x2_chb, 0, 6, 1, 1)
        self.sn3_fa16x0_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_5)
        self.sn3_fa16x0_chb.setChecked(True)
        self.sn3_fa16x0_chb.setObjectName("sn3_fa16x0_chb")
        self.gridLayout_8.addWidget(self.sn3_fa16x0_chb, 0, 2, 1, 1)
        self.sn3_fa12x0_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_5)
        self.sn3_fa12x0_chb.setObjectName("sn3_fa12x0_chb")
        self.gridLayout_8.addWidget(self.sn3_fa12x0_chb, 0, 0, 1, 1)
        self.sn3_fa20x3_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_5)
        self.sn3_fa20x3_chb.setChecked(True)
        self.sn3_fa20x3_chb.setObjectName("sn3_fa20x3_chb")
        self.gridLayout_8.addWidget(self.sn3_fa20x3_chb, 1, 1, 1, 1)
        self.lipidgen_tab1_sn3_verticalLayout.addLayout(self.gridLayout_8)
        self.horizontalLayout_33 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_33.setObjectName("horizontalLayout_33")
        self.label_61 = QtWidgets.QLabel(self.verticalLayoutWidget_5)
        self.label_61.setObjectName("label_61")
        self.horizontalLayout_33.addWidget(self.label_61)
        self.sn3_otherfa_lb = QtWidgets.QLineEdit(self.verticalLayoutWidget_5)
        self.sn3_otherfa_lb.setObjectName("sn3_otherfa_lb")
        self.horizontalLayout_33.addWidget(self.sn3_otherfa_lb)
        self.lipidgen_tab1_sn3_verticalLayout.addLayout(self.horizontalLayout_33)
        self.line_2 = QtWidgets.QFrame(self.lipidgen_tab1)
        self.line_2.setGeometry(QtCore.QRect(20, 280 * scale, 861 * scale, 16 * scale))
        self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName("line_2")
        self.line_3 = QtWidgets.QFrame(self.lipidgen_tab1)
        self.line_3.setGeometry(QtCore.QRect(20, 170 * scale, 861 * scale, 16 * scale))
        self.line_3.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_3.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_3.setObjectName("line_3")
        self.line_4 = QtWidgets.QFrame(self.lipidgen_tab1)
        self.line_4.setGeometry(QtCore.QRect(20, 20, 861 * scale, 16 * scale))
        self.line_4.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_4.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_4.setObjectName("line_4")
        self.verticalLayoutWidget_6 = QtWidgets.QWidget(self.lipidgen_tab1)
        self.verticalLayoutWidget_6.setGeometry(
            QtCore.QRect(20, 400 * scale, 859 * scale, 100 * scale)
        )
        self.verticalLayoutWidget_6.setObjectName("verticalLayoutWidget_6")
        self.lipidgen_tab1_sn4_verticalLayout = QtWidgets.QVBoxLayout(
            self.verticalLayoutWidget_6
        )
        self.lipidgen_tab1_sn4_verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.lipidgen_tab1_sn4_verticalLayout.setObjectName(
            "lipidgen_tab1_sn4_verticalLayout"
        )
        self.horizontalLayout_34 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_34.setObjectName("horizontalLayout_34")
        self.lipidgen_tab1_sn4_lb = QtWidgets.QLabel(self.verticalLayoutWidget_6)
        self.lipidgen_tab1_sn4_lb.setObjectName("lipidgen_tab1_sn4_lb")
        self.horizontalLayout_34.addWidget(self.lipidgen_tab1_sn4_lb)
        spacerItem34 = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_34.addItem(spacerItem34)
        self.lipidgen_tab1_sn3allfa_chb_2 = QtWidgets.QCheckBox(
            self.verticalLayoutWidget_6
        )
        self.lipidgen_tab1_sn3allfa_chb_2.setObjectName("lipidgen_tab1_sn3allfa_chb_2")
        self.horizontalLayout_34.addWidget(self.lipidgen_tab1_sn3allfa_chb_2)
        self.lipidgen_tab1_sn4_verticalLayout.addLayout(self.horizontalLayout_34)
        self.gridLayout_10 = QtWidgets.QGridLayout()
        self.gridLayout_10.setObjectName("gridLayout_10")
        self.sn4_fa22x6_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_6)
        self.sn4_fa22x6_chb.setChecked(True)
        self.sn4_fa22x6_chb.setObjectName("sn4_fa22x6_chb")
        self.gridLayout_10.addWidget(self.sn4_fa22x6_chb, 1, 6, 1, 1)
        self.sn4_fa22x5_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_6)
        self.sn4_fa22x5_chb.setChecked(True)
        self.sn4_fa22x5_chb.setObjectName("sn4_fa22x5_chb")
        self.gridLayout_10.addWidget(self.sn4_fa22x5_chb, 1, 5, 1, 1)
        self.sn4_fa22x4_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_6)
        self.sn4_fa22x4_chb.setChecked(True)
        self.sn4_fa22x4_chb.setObjectName("sn4_fa22x4_chb")
        self.gridLayout_10.addWidget(self.sn4_fa22x4_chb, 1, 4, 1, 1)
        self.sn4_fa20x5_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_6)
        self.sn4_fa20x5_chb.setChecked(True)
        self.sn4_fa20x5_chb.setObjectName("sn4_fa20x5_chb")
        self.gridLayout_10.addWidget(self.sn4_fa20x5_chb, 1, 3, 1, 1)
        self.sn4_fa20x4_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_6)
        self.sn4_fa20x4_chb.setChecked(True)
        self.sn4_fa20x4_chb.setObjectName("sn4_fa20x4_chb")
        self.gridLayout_10.addWidget(self.sn4_fa20x4_chb, 1, 2, 1, 1)
        self.sn4_fa20x3_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_6)
        self.sn4_fa20x3_chb.setChecked(True)
        self.sn4_fa20x3_chb.setObjectName("sn4_fa20x3_chb")
        self.gridLayout_10.addWidget(self.sn4_fa20x3_chb, 1, 1, 1, 1)
        self.sn4_fa18x3_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_6)
        self.sn4_fa18x3_chb.setChecked(True)
        self.sn4_fa18x3_chb.setObjectName("sn4_fa18x3_chb")
        self.gridLayout_10.addWidget(self.sn4_fa18x3_chb, 1, 0, 1, 1)
        self.sn4_fa18x2_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_6)
        self.sn4_fa18x2_chb.setChecked(True)
        self.sn4_fa18x2_chb.setObjectName("sn4_fa18x2_chb")
        self.gridLayout_10.addWidget(self.sn4_fa18x2_chb, 0, 6, 1, 1)
        self.sn4_fa18x1_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_6)
        self.sn4_fa18x1_chb.setChecked(True)
        self.sn4_fa18x1_chb.setObjectName("sn4_fa18x1_chb")
        self.gridLayout_10.addWidget(self.sn4_fa18x1_chb, 0, 5, 1, 1)
        self.sn4_fa18x0_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_6)
        self.sn4_fa18x0_chb.setChecked(True)
        self.sn4_fa18x0_chb.setObjectName("sn4_fa18x0_chb")
        self.gridLayout_10.addWidget(self.sn4_fa18x0_chb, 0, 4, 1, 1)
        self.sn4_fa16x1_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_6)
        self.sn4_fa16x1_chb.setChecked(False)
        self.sn4_fa16x1_chb.setObjectName("sn4_fa16x1_chb")
        self.gridLayout_10.addWidget(self.sn4_fa16x1_chb, 0, 3, 1, 1)
        self.sn4_fa16x0_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_6)
        self.sn4_fa16x0_chb.setChecked(True)
        self.sn4_fa16x0_chb.setObjectName("sn4_fa16x0_chb")
        self.gridLayout_10.addWidget(self.sn4_fa16x0_chb, 0, 2, 1, 1)
        self.sn4_fa14x0_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_6)
        self.sn4_fa14x0_chb.setObjectName("sn4_fa14x0_chb")
        self.gridLayout_10.addWidget(self.sn4_fa14x0_chb, 0, 1, 1, 1)
        self.sn4_fa12x0_chb = QtWidgets.QCheckBox(self.verticalLayoutWidget_6)
        self.sn4_fa12x0_chb.setObjectName("sn4_fa12x0_chb")
        self.gridLayout_10.addWidget(self.sn4_fa12x0_chb, 0, 0, 1, 1)
        self.lipidgen_tab1_sn4_verticalLayout.addLayout(self.gridLayout_10)
        self.horizontalLayout_35 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_35.setObjectName("horizontalLayout_35")
        self.label_66 = QtWidgets.QLabel(self.verticalLayoutWidget_6)
        self.label_66.setObjectName("label_66")
        self.horizontalLayout_35.addWidget(self.label_66)
        self.sn4_otherfa_lb = QtWidgets.QLineEdit(self.verticalLayoutWidget_6)
        self.sn4_otherfa_lb.setObjectName("sn4_otherfa_lb")
        self.horizontalLayout_35.addWidget(self.sn4_otherfa_lb)
        self.lipidgen_tab1_sn4_verticalLayout.addLayout(self.horizontalLayout_35)
        self.line_7 = QtWidgets.QFrame(self.lipidgen_tab1)
        self.line_7.setGeometry(QtCore.QRect(20, 390 * scale, 861 * scale, 16 * scale))
        self.line_7.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_7.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_7.setObjectName("line_7")
        self.line_11 = QtWidgets.QFrame(self.lipidgen_tab1)
        self.line_11.setGeometry(QtCore.QRect(20, 500 * scale, 861 * scale, 16 * scale))
        self.line_11.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_11.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_11.setObjectName("line_11")
        self.horizontalLayoutWidget_6 = QtWidgets.QWidget(self.lipidgen_tab1)
        self.horizontalLayoutWidget_6.setGeometry(
            QtCore.QRect(20, 540 * scale, 781 * scale, 31 * scale)
        )
        self.horizontalLayoutWidget_6.setObjectName("horizontalLayoutWidget_6")
        self.horizontalLayout_40 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_6)
        self.horizontalLayout_40.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_40.setObjectName("horizontalLayout_40")
        self.lipidgen_tab1_genlist_pb = QtWidgets.QPushButton(
            self.horizontalLayoutWidget_6
        )
        self.lipidgen_tab1_genlist_pb.setMinimumSize(QtCore.QSize(120, 0))
        self.lipidgen_tab1_genlist_pb.setObjectName("lipidgen_tab1_genlist_pb")
        self.horizontalLayout_40.addWidget(self.lipidgen_tab1_genlist_pb)
        self.lipidgen_tab1_genlist_le = QtWidgets.QLineEdit(
            self.horizontalLayoutWidget_6
        )
        self.lipidgen_tab1_genlist_le.setObjectName("lipidgen_tab1_genlist_le")
        self.horizontalLayout_40.addWidget(self.lipidgen_tab1_genlist_le)
        self.lipidgen_tabframe.addTab(self.lipidgen_tab1, "")
        self.lipidgen_tab2 = QtWidgets.QWidget()
        self.lipidgen_tab2.setObjectName("lipidgen_tab2")
        self.horizontalLayoutWidget_10 = QtWidgets.QWidget(self.lipidgen_tab2)
        self.horizontalLayoutWidget_10.setGeometry(
            QtCore.QRect(10, 10, 881 * scale, 31 * scale)
        )
        self.horizontalLayoutWidget_10.setObjectName("horizontalLayoutWidget_10")
        self.horizontalLayout_9 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_10)
        self.horizontalLayout_9.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_9.setObjectName("horizontalLayout_9")
        self.label_60 = QtWidgets.QLabel(self.horizontalLayoutWidget_10)
        self.label_60.setObjectName("label_60")
        self.horizontalLayout_9.addWidget(self.label_60)
        self.lipidgen_tab2_loadlist_le = QtWidgets.QLineEdit(
            self.horizontalLayoutWidget_10
        )
        self.lipidgen_tab2_loadlist_le.setObjectName("lipidgen_tab2_loadlist_le")
        self.horizontalLayout_9.addWidget(self.lipidgen_tab2_loadlist_le)
        self.lipidgen_tab2_loadlist_pb = QtWidgets.QPushButton(
            self.horizontalLayoutWidget_10
        )
        self.lipidgen_tab2_loadlist_pb.setObjectName("lipidgen_tab2_loadlist_pb")
        self.horizontalLayout_9.addWidget(self.lipidgen_tab2_loadlist_pb)
        self.lipidgen_tab2_loadlist_lb = QtWidgets.QLabel(self.lipidgen_tab2)
        self.lipidgen_tab2_loadlist_lb.setGeometry(
            QtCore.QRect(10, 50 * scale, 881 * scale, 16 * scale)
        )
        self.lipidgen_tab2_loadlist_lb.setObjectName("lipidgen_tab2_loadlist_lb")
        self.horizontalLayoutWidget_11 = QtWidgets.QWidget(self.lipidgen_tab2)
        self.horizontalLayoutWidget_11.setGeometry(
            QtCore.QRect(10, 70 * scale, 881 * scale, 31 * scale)
        )
        self.horizontalLayoutWidget_11.setObjectName("horizontalLayoutWidget_11")
        self.horizontalLayout_41 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_11)
        self.horizontalLayout_41.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_41.setObjectName("horizontalLayout_41")
        self.lipidgen_tab2_savemasterlist_lb = QtWidgets.QLabel(
            self.horizontalLayoutWidget_11
        )
        self.lipidgen_tab2_savemasterlist_lb.setObjectName(
            "lipidgen_tab2_savemasterlist_lb"
        )
        self.horizontalLayout_41.addWidget(self.lipidgen_tab2_savemasterlist_lb)
        self.lipidgen_tab2_savemasterlist_le = QtWidgets.QLineEdit(
            self.horizontalLayoutWidget_11
        )
        self.lipidgen_tab2_savemasterlist_le.setObjectName(
            "lipidgen_tab2_savemasterlist_le"
        )
        self.horizontalLayout_41.addWidget(self.lipidgen_tab2_savemasterlist_le)
        self.lipidgen_tab2_savemasterlist_pb = QtWidgets.QPushButton(
            self.horizontalLayoutWidget_11
        )
        self.lipidgen_tab2_savemasterlist_pb.setObjectName(
            "lipidgen_tab2_savemasterlist_pb"
        )
        self.horizontalLayout_41.addWidget(self.lipidgen_tab2_savemasterlist_pb)
        self.horizontalLayoutWidget_8 = QtWidgets.QWidget(self.lipidgen_tab2)
        self.horizontalLayoutWidget_8.setGeometry(
            QtCore.QRect(10, 110 * scale, 881 * scale, 231 * scale)
        )
        self.horizontalLayoutWidget_8.setObjectName("horizontalLayoutWidget_8")
        self.horizontalLayout_42 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_8)
        self.horizontalLayout_42.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_42.setObjectName("horizontalLayout_42")
        self.lipidgen_tab1_genmasterlist_pb = QtWidgets.QPushButton(
            self.horizontalLayoutWidget_8
        )
        self.lipidgen_tab1_genmasterlist_pb.setMinimumSize(QtCore.QSize(120, 0))
        self.lipidgen_tab1_genmasterlist_pb.setObjectName(
            "lipidgen_tab1_genmasterlist_pb"
        )
        self.horizontalLayout_42.addWidget(self.lipidgen_tab1_genmasterlist_pb)
        self.lipidgen_tab1_genmasterlist_pte = QtWidgets.QPlainTextEdit(
            self.horizontalLayoutWidget_8
        )
        self.lipidgen_tab1_genmasterlist_pte.setObjectName(
            "lipidgen_tab1_genmasterlist_pte"
        )
        self.horizontalLayout_42.addWidget(self.lipidgen_tab1_genmasterlist_pte)
        self.lipidgen_tabframe.addTab(self.lipidgen_tab2, "")
        self.tabframe.addTab(self.generator_tab, "")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.verticalLayoutWidget = QtWidgets.QWidget(self.tab)
        self.verticalLayoutWidget.setGeometry(
            QtCore.QRect(210 * scale, 10, 701 * scale, 201 * scale)
        )
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.version_lb = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.version_lb.setObjectName("version_lb")
        self.verticalLayout_4.addWidget(self.version_lb)
        self.label_26 = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.label_26.setObjectName("label_26")
        self.verticalLayout_4.addWidget(self.label_26)
        self.label_43 = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.label_43.setObjectName("label_43")
        self.verticalLayout_4.addWidget(self.label_43)
        self.formLayout = QtWidgets.QFormLayout()
        self.formLayout.setObjectName("formLayout")
        self.label_42 = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.label_42.setObjectName("label_42")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_42)
        self.label_44 = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.label_44.setObjectName("label_44")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_44)
        self.link_gplv2_lb = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.link_gplv2_lb.setObjectName("link_gplv2_lb")
        self.formLayout.setWidget(
            0, QtWidgets.QFormLayout.FieldRole, self.link_gplv2_lb
        )
        self.label_45 = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.label_45.setTextFormat(QtCore.Qt.RichText)
        self.label_45.setObjectName("label_45")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.label_45)
        self.verticalLayout_4.addLayout(self.formLayout)
        self.link_source_lb = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.link_source_lb.setObjectName("link_source_lb")
        self.verticalLayout_4.addWidget(self.link_source_lb)
        self.link_tutorial_lb = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.link_tutorial_lb.setObjectName("link_tutorial_lb")
        self.verticalLayout_4.addWidget(self.link_tutorial_lb)
        self.link_paper_lb = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.link_paper_lb.setWordWrap(True)
        self.link_paper_lb.setObjectName("link_paper_lb")
        self.verticalLayout_4.addWidget(self.link_paper_lb)
        self.link_otherprojects_lb = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.link_otherprojects_lb.setObjectName("link_otherprojects_lb")
        self.verticalLayout_4.addWidget(self.link_otherprojects_lb)
        self.line_8 = QtWidgets.QFrame(self.tab)
        self.line_8.setGeometry(QtCore.QRect(10, 210 * scale, 901 * scale, 16 * scale))
        self.line_8.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_8.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_8.setObjectName("line_8")
        self.label_30 = QtWidgets.QLabel(self.tab)
        self.label_30.setGeometry(
            QtCore.QRect(20, 228 * scale, 881 * scale, 16 * scale)
        )
        self.label_30.setObjectName("label_30")
        self.logo_pb = QtWidgets.QPushButton(self.tab)
        self.logo_pb.setGeometry(QtCore.QRect(40, 40, 128 * scale, 128 * scale))
        self.logo_pb.setMinimumSize(QtCore.QSize(128, 128))
        self.logo_pb.setMaximumSize(QtCore.QSize(128, 128))
        self.logo_pb.setText("")
        icon1 = QtGui.QIcon()
        icon1.addPixmap(
            QtGui.QPixmap(":/LipidHunter.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off
        )
        self.logo_pb.setIcon(icon1)
        self.logo_pb.setIconSize(QtCore.QSize(128, 128))
        self.logo_pb.setObjectName("logo_pb")
        self.formLayoutWidget_3 = QtWidgets.QWidget(self.tab)
        self.formLayoutWidget_3.setGeometry(
            QtCore.QRect(20, 250 * scale, 881 * scale, 261 * scale)
        )
        self.formLayoutWidget_3.setObjectName("formLayoutWidget_3")
        self.formLayout_3 = QtWidgets.QFormLayout(self.formLayoutWidget_3)
        self.formLayout_3.setFieldGrowthPolicy(
            QtWidgets.QFormLayout.AllNonFixedFieldsGrow
        )
        self.formLayout_3.setFormAlignment(
            QtCore.Qt.AlignLeading | QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop
        )
        self.formLayout_3.setContentsMargins(0, 0, 0, 0)
        self.formLayout_3.setObjectName("formLayout_3")
        self.label_31 = QtWidgets.QLabel(self.formLayoutWidget_3)
        self.label_31.setObjectName("label_31")
        self.formLayout_3.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_31)
        self.link_matplotlib_lb = QtWidgets.QLabel(self.formLayoutWidget_3)
        self.link_matplotlib_lb.setObjectName("link_matplotlib_lb")
        self.formLayout_3.setWidget(
            0, QtWidgets.QFormLayout.FieldRole, self.link_matplotlib_lb
        )
        self.label_36 = QtWidgets.QLabel(self.formLayoutWidget_3)
        self.label_36.setObjectName("label_36")
        self.formLayout_3.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.label_36)
        self.link_numpy_lb = QtWidgets.QLabel(self.formLayoutWidget_3)
        self.link_numpy_lb.setMinimumSize(QtCore.QSize(0, 15))
        self.link_numpy_lb.setMaximumSize(QtCore.QSize(800, 16777215))
        self.link_numpy_lb.setTextFormat(QtCore.Qt.RichText)
        self.link_numpy_lb.setWordWrap(True)
        self.link_numpy_lb.setObjectName("link_numpy_lb")
        self.formLayout_3.setWidget(
            2, QtWidgets.QFormLayout.FieldRole, self.link_numpy_lb
        )
        self.label_32 = QtWidgets.QLabel(self.formLayoutWidget_3)
        self.label_32.setObjectName("label_32")
        self.formLayout_3.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.label_32)
        self.link_pandas_lb = QtWidgets.QLabel(self.formLayoutWidget_3)
        self.link_pandas_lb.setObjectName("link_pandas_lb")
        self.formLayout_3.setWidget(
            3, QtWidgets.QFormLayout.FieldRole, self.link_pandas_lb
        )
        self.label_33 = QtWidgets.QLabel(self.formLayoutWidget_3)
        self.label_33.setObjectName("label_33")
        self.formLayout_3.setWidget(4, QtWidgets.QFormLayout.LabelRole, self.label_33)
        self.link_pymzml_lb = QtWidgets.QLabel(self.formLayoutWidget_3)
        self.link_pymzml_lb.setMinimumSize(QtCore.QSize(0, 15))
        self.link_pymzml_lb.setMaximumSize(QtCore.QSize(800, 16777215))
        self.link_pymzml_lb.setTextFormat(QtCore.Qt.RichText)
        self.link_pymzml_lb.setWordWrap(True)
        self.link_pymzml_lb.setObjectName("link_pymzml_lb")
        self.formLayout_3.setWidget(
            4, QtWidgets.QFormLayout.FieldRole, self.link_pymzml_lb
        )
        self.label_34 = QtWidgets.QLabel(self.formLayoutWidget_3)
        self.label_34.setObjectName("label_34")
        self.formLayout_3.setWidget(5, QtWidgets.QFormLayout.LabelRole, self.label_34)
        self.link_pyside_lb = QtWidgets.QLabel(self.formLayoutWidget_3)
        self.link_pyside_lb.setObjectName("link_pyside_lb")
        self.formLayout_3.setWidget(
            5, QtWidgets.QFormLayout.FieldRole, self.link_pyside_lb
        )
        self.label_37 = QtWidgets.QLabel(self.formLayoutWidget_3)
        self.label_37.setObjectName("label_37")
        self.formLayout_3.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_37)
        self.label_97 = QtWidgets.QLabel(self.formLayoutWidget_3)
        self.label_97.setObjectName("label_97")
        self.formLayout_3.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.label_97)
        self.line_10 = QtWidgets.QFrame(self.tab)
        self.line_10.setGeometry(QtCore.QRect(10, 510 * scale, 901 * scale, 16 * scale))
        self.line_10.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_10.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_10.setObjectName("line_10")
        self.label_46 = QtWidgets.QLabel(self.tab)
        self.label_46.setGeometry(
            QtCore.QRect(20 * scale, 520 * scale, 881 * scale, 16 * scale)
        )
        self.label_46.setObjectName("label_46")
        self.formLayoutWidget_2 = QtWidgets.QWidget(self.tab)
        self.formLayoutWidget_2.setGeometry(
            QtCore.QRect(20 * scale, 540 * scale, 881 * scale, 71 * scale)
        )
        self.formLayoutWidget_2.setObjectName("formLayoutWidget_2")
        self.formLayout_2 = QtWidgets.QFormLayout(self.formLayoutWidget_2)
        self.formLayout_2.setFormAlignment(
            QtCore.Qt.AlignLeading | QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop
        )
        self.formLayout_2.setContentsMargins(0, 0, 0, 0)
        self.formLayout_2.setObjectName("formLayout_2")
        self.label_47 = QtWidgets.QLabel(self.formLayoutWidget_2)
        self.label_47.setObjectName("label_47")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_47)
        self.label_49 = QtWidgets.QLabel(self.formLayoutWidget_2)
        self.label_49.setAlignment(
            QtCore.Qt.AlignRight | QtCore.Qt.AlignTrailing | QtCore.Qt.AlignVCenter
        )
        self.label_49.setObjectName("label_49")
        self.formLayout_2.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_49)
        self.link_bmbf_lb = QtWidgets.QLabel(self.formLayoutWidget_2)
        self.link_bmbf_lb.setObjectName("link_bmbf_lb")
        self.formLayout_2.setWidget(
            0, QtWidgets.QFormLayout.FieldRole, self.link_bmbf_lb
        )
        self.link_emed_lb = QtWidgets.QLabel(self.formLayoutWidget_2)
        self.link_emed_lb.setObjectName("link_emed_lb")
        self.formLayout_2.setWidget(
            1, QtWidgets.QFormLayout.FieldRole, self.link_emed_lb
        )
        self.link_sysmedos_lb = QtWidgets.QLabel(self.formLayoutWidget_2)
        self.link_sysmedos_lb.setObjectName("link_sysmedos_lb")
        self.formLayout_2.setWidget(
            2, QtWidgets.QFormLayout.FieldRole, self.link_sysmedos_lb
        )
        self.label_58 = QtWidgets.QLabel(self.formLayoutWidget_2)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_58.sizePolicy().hasHeightForWidth())
        self.label_58.setSizePolicy(sizePolicy)
        self.label_58.setAlignment(
            QtCore.Qt.AlignRight | QtCore.Qt.AlignTrailing | QtCore.Qt.AlignVCenter
        )
        self.label_58.setObjectName("label_58")
        self.formLayout_2.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.label_58)
        self.tabframe.addTab(self.tab, "")
        self.link_uni_lb = QtWidgets.QLabel(self.centralwidget)
        self.link_uni_lb.setGeometry(
            QtCore.QRect(20 * scale, 650 * scale, 901 * scale, 16 * scale)
        )
        self.link_uni_lb.setContextMenuPolicy(QtCore.Qt.NoContextMenu)
        self.link_uni_lb.setObjectName("link_uni_lb")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 960 * scale, 21 * scale))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        self.tabframe.setCurrentIndex(0)
        self.tab_a_lipidclass_cmb.setCurrentIndex(0)
        self.runhunter_tabframe.setCurrentIndex(0)
        self.tab_c_lipidclass_cmb.setCurrentIndex(0)
        self.lipidgen_tabframe.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(
            QtWidgets.QApplication.translate(
                "MainWindow", "LipidHunter 2 (RC2) | University of Leipzig ", None
            )
        )
        self.label_50.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Score weight factor settings:", None
            )
        )
        self.tab_a_mzml_pb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Open", None)
        )
        self.label_3.setText(
            QtWidgets.QApplication.translate("MainWindow", "Select lipid class:", None)
        )
        self.label_4.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Scan time range (min):", None
            )
        )
        self.label_9.setText(QtWidgets.QApplication.translate("MainWindow", "to", None))
        self.label_10.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                '<html><head/><body><p><span style=" font-style:italic;">m/z</span>&nbsp; range:</p></body></html>',
                None,
            )
        )
        self.label_11.setText(
            QtWidgets.QApplication.translate("MainWindow", "to", None)
        )
        self.label_6.setText(
            QtWidgets.QApplication.translate("MainWindow", "Input .mzML file:", None)
        )
        self.tab_a_loadfalist_pb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Open", None)
        )
        self.label_7.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Save report to folder:", None
            )
        )
        self.tab_a_loadscorecfg_pb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Open", None)
        )
        self.tab_a_sumxlsxpath_pb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Save as", None)
        )
        self.tab_a_saveimgfolder_pb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Save as", None)
        )
        self.label_8.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Export summary .xlsx table:", None
            )
        )
        self.label_16.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "MS/MS level tolerance:", None
            )
        )
        self.tab_a_ms2ppm_spb.setSuffix(
            QtWidgets.QApplication.translate("MainWindow", " ppm", None)
        )
        self.tab_a_ms2ppm_spb.setPrefix(
            QtWidgets.QApplication.translate("MainWindow", "+/- ", None)
        )
        self.tab_d_ms2threshold_spb_2.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "MS/MS level threshold (absolute intensity):", None
            )
        )
        self.tab_a_lipidclass_cmb.setItemText(
            0,
            QtWidgets.QApplication.translate(
                "MainWindow", "--- Please select ---", None
            ),
        )
        self.tab_a_lipidclass_cmb.setItemText(
            1,
            QtWidgets.QApplication.translate(
                "MainWindow", "Phosphatidic acid (PA) [M-H]-", None
            ),
        )
        self.tab_a_lipidclass_cmb.setItemText(
            2,
            QtWidgets.QApplication.translate(
                "MainWindow", "Phosphatidylcholine (PC) [M+HCOO]-", None
            ),
        )
        self.tab_a_lipidclass_cmb.setItemText(
            3,
            QtWidgets.QApplication.translate(
                "MainWindow", "Phosphatidylcholine (PC) [M+CH3COO]-", None
            ),
        )
        self.tab_a_lipidclass_cmb.setItemText(
            4,
            QtWidgets.QApplication.translate(
                "MainWindow", "Phosphatidylethanolamine (PE) [M-H]-", None
            ),
        )
        self.tab_a_lipidclass_cmb.setItemText(
            5,
            QtWidgets.QApplication.translate(
                "MainWindow", "Phosphatidylglycerol (PG) [M-H]-", None
            ),
        )
        self.tab_a_lipidclass_cmb.setItemText(
            6,
            QtWidgets.QApplication.translate(
                "MainWindow", "Phosphatidylinositol (PI) [M-H]-", None
            ),
        )
        self.tab_a_lipidclass_cmb.setItemText(
            7,
            QtWidgets.QApplication.translate(
                "MainWindow", "Phosphatidylserine (PS) [M-H]-", None
            ),
        )
        self.tab_a_lipidclass_cmb.setItemText(
            8,
            QtWidgets.QApplication.translate(
                "MainWindow", "Triacylglycerol (TG) [M+NH4]+", None
            ),
        )
        self.tab_a_lipidclass_cmb.setItemText(
            9,
            QtWidgets.QApplication.translate(
                "MainWindow", "Triacylglycerol (TG) [M+H]+", None
            ),
        )
        self.tab_a_lipidclass_cmb.setItemText(
            10,
            QtWidgets.QApplication.translate(
                "MainWindow", "Triacylglycerol (TG) [M+Na]+", None
            ),
        )
        self.tab_a_lipidclass_cmb.setItemText(
            11,
            QtWidgets.QApplication.translate(
                "MainWindow", "Diacylglycerol (DG) [M+NH4]+", None
            ),
        )
        # self.tab_a_lipidclass_cmb.setItemText(12, QtWidgets.QApplication.translate("MainWindow", "Diacylglycerol (DG) [M+H]+", None, ))
        self.tab_a_lipidclass_cmb.setItemText(
            13 - 1,
            QtWidgets.QApplication.translate(
                "MainWindow", "Lyso PA (LPA) [M-H]-", None
            ),
        )
        self.tab_a_lipidclass_cmb.setItemText(
            14 - 1,
            QtWidgets.QApplication.translate(
                "MainWindow", "Lyso PC (LPC) [M+HCOO]-", None
            ),
        )
        self.tab_a_lipidclass_cmb.setItemText(
            15 - 1,
            QtWidgets.QApplication.translate(
                "MainWindow", "Lyso PC (LPC) [M+CH3COO]-", None
            ),
        )
        self.tab_a_lipidclass_cmb.setItemText(
            16 - 1,
            QtWidgets.QApplication.translate(
                "MainWindow", "Lyso PE (LPE) [M-H]-", None
            ),
        )
        self.tab_a_lipidclass_cmb.setItemText(
            17 - 1,
            QtWidgets.QApplication.translate(
                "MainWindow", "Lyso PG (LPG) [M-H]-", None
            ),
        )
        self.tab_a_lipidclass_cmb.setItemText(
            18 - 1,
            QtWidgets.QApplication.translate(
                "MainWindow", "Lyso PI (LPI) [M-H]-", None
            ),
        )
        self.tab_a_lipidclass_cmb.setItemText(
            19 - 1,
            QtWidgets.QApplication.translate(
                "MainWindow", "Lyso PS (LPS) [M-H]-", None
            ),
        )
        self.label_5.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA white list:", None)
        )
        self.label_19.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                "<html><head/><body><p>Isotope score &ge;</p></body></html>",
                None,
            )
        )
        self.label_20.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                "<html><head/><body><p>Rank score &ge;</p></body></html>",
                None,
            )
        )
        self.label_38.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Minimum relative intensity for scoring:", None
            )
        )
        self.tab_a_ms2infoth_dspb.setSuffix(
            QtWidgets.QApplication.translate("MainWindow", " %", None)
        )
        self.label_15.setText(
            QtWidgets.QApplication.translate("MainWindow", "MS level tolerance:", None)
        )
        self.tab_a_msppm_spb.setSuffix(
            QtWidgets.QApplication.translate("MainWindow", " ppm", None)
        )
        self.tab_a_msppm_spb.setPrefix(
            QtWidgets.QApplication.translate("MainWindow", "+/- ", None)
        )
        self.label_13.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "MS level threshold (absolute intensity):", None
            )
        )
        self.tab_a_msmax_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Set MS max inetnsity", None)
        )
        self.label_25.setText(
            QtWidgets.QApplication.translate("MainWindow", "Precursor window:", None)
        )
        self.tab_a_prwindow_spb.setPrefix(
            QtWidgets.QApplication.translate("MainWindow", "+/- ", None)
        )
        self.tab_a_prwindow_spb.setSuffix(
            QtWidgets.QApplication.translate("MainWindow", " m/z", None)
        )
        self.x.setText(QtWidgets.QApplication.translate("MainWindow", "DDA Top:", None))
        self.tab_a_runhunter_pb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Hunt for lipids !", None)
        )
        self.runhunter_tabframe.setTabText(
            self.runhunter_tabframe.indexOf(self.singlerun_tab),
            QtWidgets.QApplication.translate(
                "MainWindow", "Run LipidHunter identification directly", None
            ),
        )
        self.label.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Save configuration file as:", None
            )
        )
        self.tab_a_cfgpath_pb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Save As", None)
        )
        self.tab_a_gencfg_pb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Generate configuration file for batch mode >>>", None
            )
        )
        self.runhunter_tabframe.setTabText(
            self.runhunter_tabframe.indexOf(self.savecfg_tab),
            QtWidgets.QApplication.translate(
                "MainWindow", "Save parameters for Batch Mode Hunter", None
            ),
        )
        self.label_23.setText(
            QtWidgets.QApplication.translate("MainWindow", "Instrument vendor:", None)
        )
        self.vendor_cmb.setItemText(
            0,
            QtWidgets.QApplication.translate(
                "MainWindow", "--- Please select ---", None
            ),
        )
        self.vendor_cmb.setItemText(
            1, QtWidgets.QApplication.translate("MainWindow", "Agilent", None)
        )
        self.vendor_cmb.setItemText(
            2, QtWidgets.QApplication.translate("MainWindow", "SCIEX", None)
        )
        self.vendor_cmb.setItemText(
            3, QtWidgets.QApplication.translate("MainWindow", "Thermo", None)
        )
        self.vendor_cmb.setItemText(
            4, QtWidgets.QApplication.translate("MainWindow", "Waters", None)
        )
        self.label_24.setText(
            QtWidgets.QApplication.translate("MainWindow", "Mode:", None)
        )
        self.mode_lcms_rb.setText(
            QtWidgets.QApplication.translate("MainWindow", "LC-MS/MS", None)
        )
        self.mode_static_rb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Shotgun", None)
        )
        self.tabframe.setTabText(
            self.tabframe.indexOf(self.hunter_tab),
            QtWidgets.QApplication.translate("MainWindow", "Hunt for Lipids", None),
        )
        self.label_98.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Load batch configurations files (.txt):", None
            )
        )
        self.tab_b_1_lb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Input batch configuration file:", None
            )
        )
        self.tab_b_addcfg_pb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Add a single file", None)
        )
        self.tab_b_2_lb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Add all batch configurations in a folder:", None
            )
        )
        self.tab_b_addcfgfolder_pb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Select a folder", None)
        )
        self.tab_b_clearall_pb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Clear all", None)
        )
        self.tab_b_runbatch_pb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Run batch mode identification >>>", None
            )
        )
        self.tab_b_maxsubcore_lb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Max cores for each task:", None
            )
        )
        self.tab_b_maxsubram_lb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Max RAM for each task:", None
            )
        )
        self.tab_b_maxsubram_spb.setSuffix(
            QtWidgets.QApplication.translate("MainWindow", " GB", None)
        )
        self.tab_b_mutlimode_cmb.setItemText(
            0,
            QtWidgets.QApplication.translate(
                "MainWindow", "Process as a sequence", None
            ),
        )
        self.tab_b_mutlimode_cmb.setItemText(
            1,
            QtWidgets.QApplication.translate(
                "MainWindow", "Process multiple files in parallel", None
            ),
        )
        self.tab_b_maxbatch_lb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Max number of files in parallel:", None
            )
        )
        self.label_99.setText(
            QtWidgets.QApplication.translate("MainWindow", "Batch mode settings:", None)
        )
        self.tabframe.setTabText(
            self.tabframe.indexOf(self.batch_tab),
            QtWidgets.QApplication.translate("MainWindow", "Batch Mode Hunter", None),
        )
        self.label_93.setText(
            QtWidgets.QApplication.translate("MainWindow", "Image format:", None)
        )
        self.tab_c_imagetype_cmb.setItemText(
            0, QtWidgets.QApplication.translate("MainWindow", ".png", None)
        )
        self.tab_c_imagetype_cmb.setItemText(
            1, QtWidgets.QApplication.translate("MainWindow", ".svg", None)
        )
        self.label_94.setText(
            QtWidgets.QApplication.translate("MainWindow", "dpi:", None)
        )
        self.tab_c_scorecfgtg_lb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                '<html><head/><body><p>&quot;W<span style=" font-size:11pt; vertical-align:sub;">frag</span>&quot;  for TG:</p></body></html>',
                None,
            )
        )
        self.tab_c_scorecfgdg_lb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                '<html><head/><body><p>&quot;W<span style=" font-size:11pt; vertical-align:sub;">frag</span>&quot;  for DG:</p></body></html>',
                None,
            )
        )
        self.tab_c_isotopescoremode_lb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Isotope Score mode:", None)
        )
        self.label_14.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Lipid specific FRAG & NL list:", None
            )
        )
        self.tab_c_scorecfgdg_pb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Open", None)
        )
        self.tab_c_scorecfgpl_lb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                '<html><head/><body><p>&quot;W<span style=" font-size:11pt; vertical-align:sub;">frag</span>&quot; for PL:</p></body></html>',
                None,
            )
        )
        self.tab_c_scorecfgpl_pb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Open", None)
        )
        self.tab_c_scorecfgtg_pb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Open", None)
        )
        self.tab_c_hgcfg_pb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Open", None)
        )
        self.label_21.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                "<html><head/><body><p>Score weight factor settings:</p></body></html>",
                None,
            )
        )
        self.tab_c_parallization_cmb.setItemText(
            0,
            QtWidgets.QApplication.translate(
                "MainWindow", "CPU prarallization mode", None
            ),
        )
        self.tab_c_parallization_cmb.setItemText(
            1,
            QtWidgets.QApplication.translate(
                "MainWindow", "CPU +GPU mode (experimental)", None
            ),
        )
        self.label_63.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Max CPU cores to use:", None
            )
        )
        self.label_64.setText(
            QtWidgets.QApplication.translate("MainWindow", "Max RAM to use:", None)
        )
        self.tab_c_ram_spb.setSuffix(
            QtWidgets.QApplication.translate("MainWindow", " GB", None)
        )
        self.tab_c_tag_all_fa_chb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                "Favoured for TG with all three FA residues identified",
                None,
            )
        )
        self.label_12.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Default Fatty acid white list:", None
            )
        )
        self.tab_c_falistpl_pb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Open", None)
        )
        self.tab_c_savesettings_pb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Save above settings as defaullt", None
            )
        )
        self.tab_c_falisttg_lb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                '<html><head/><body><p>&quot;W<span style=" font-size:11pt; vertical-align:sub;">frag</span>&quot; for LPL:</p></body></html>',
                None,
            )
        )
        self.tab_c_falistpl_lb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Default FA list:", None)
        )
        self.tab_c_falisttg_pb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Open", None)
        )
        self.tab_c_isotopescoremode_cmb.setItemText(
            0, QtWidgets.QApplication.translate("MainWindow", "All elements", None)
        )
        self.tab_c_isotopescoremode_cmb.setItemText(
            1, QtWidgets.QApplication.translate("MainWindow", "Fast (only 13C)", None)
        )
        self.tab_c_scoremode_lb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Score mode:", None)
        )
        self.tab_c_scoremode_cmb.setItemText(
            0, QtWidgets.QApplication.translate("MainWindow", "Signal rank mode", None)
        )
        self.tab_c_scoremode_cmb.setItemText(
            1,
            QtWidgets.QApplication.translate(
                "MainWindow", "Relative Intensity mode", None
            ),
        )
        self.label_35.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Prarallel processing settings:", None
            )
        )
        self.label_54.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Output image settings:", None
            )
        )
        self.label_28.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA white list:", None)
        )
        self.tab_c_lmcalcfalist_pb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Open", None)
        )
        self.tab_c_lmexport_pb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Save as", None)
        )
        self.label_29.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                '<html><head/><body><p align="center"><span style=" font-size:10pt; font-weight:600;">Lipid Composer: Generate Lipid Master table from selected FA white list</span></p></body></html>',
                None,
            )
        )
        self.label_27.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Export Lipid Master table as (.csv): ", None
            )
        )
        self.label_48.setText(
            QtWidgets.QApplication.translate("MainWindow", "Select lipid class:", None)
        )
        self.tab_c_lipidclass_cmb.setItemText(
            0,
            QtWidgets.QApplication.translate(
                "MainWindow", "--- Please select ---", None
            ),
        )
        self.tab_c_lipidclass_cmb.setItemText(
            1,
            QtWidgets.QApplication.translate(
                "MainWindow", "Phosphatidic acid (PA) [M-H]-", None
            ),
        )
        self.tab_c_lipidclass_cmb.setItemText(
            2,
            QtWidgets.QApplication.translate(
                "MainWindow", "Phosphatidylcholine (PC) [M+HCOO]-", None
            ),
        )
        self.tab_c_lipidclass_cmb.setItemText(
            3,
            QtWidgets.QApplication.translate(
                "MainWindow", "Phosphatidylcholine (PC) [M+CH3COO]-", None
            ),
        )
        self.tab_c_lipidclass_cmb.setItemText(
            4,
            QtWidgets.QApplication.translate(
                "MainWindow", "Phosphatidylethanolamine (PE) [M-H]-", None
            ),
        )
        self.tab_c_lipidclass_cmb.setItemText(
            5,
            QtWidgets.QApplication.translate(
                "MainWindow", "Phosphatidylglycerol (PG) [M-H]-", None
            ),
        )
        self.tab_c_lipidclass_cmb.setItemText(
            6,
            QtWidgets.QApplication.translate(
                "MainWindow", "Phosphatidylinositol (PI) [M-H]-", None
            ),
        )
        self.tab_c_lipidclass_cmb.setItemText(
            7,
            QtWidgets.QApplication.translate(
                "MainWindow", "Phosphatidylserine (PS) [M-H]-", None
            ),
        )
        self.tab_c_lipidclass_cmb.setItemText(
            8,
            QtWidgets.QApplication.translate(
                "MainWindow", "Triacylglycerol (TG) [M+NH4]+", None
            ),
        )
        self.tab_c_lipidclass_cmb.setItemText(
            9,
            QtWidgets.QApplication.translate(
                "MainWindow", "Triacylglycerol (TG) [M+H]+", None
            ),
        )
        self.tab_c_lipidclass_cmb.setItemText(
            10,
            QtWidgets.QApplication.translate(
                "MainWindow", "Triacylglycerol (TG) [M+Na]+", None
            ),
        )
        self.tab_c_lipidclass_cmb.setItemText(
            11,
            QtWidgets.QApplication.translate(
                "MainWindow", "Diacylglycerol (DG) [M+NH4]+", None
            ),
        )
        # self.tab_c_lipidclass_cmb.setItemText(12, QtWidgets.QApplication.translate("MainWindow", "Diacylglycerol (DG) [M+H]+", None, ))
        self.tab_c_lipidclass_cmb.setItemText(
            13 - 1,
            QtWidgets.QApplication.translate(
                "MainWindow", "Lyso PA (LPA) [M-H]-", None
            ),
        )
        self.tab_c_lipidclass_cmb.setItemText(
            14 - 1,
            QtWidgets.QApplication.translate(
                "MainWindow", "Lyso PC (LPC) [M+HCOO]-", None
            ),
        )
        self.tab_c_lipidclass_cmb.setItemText(
            15 - 1,
            QtWidgets.QApplication.translate(
                "MainWindow", "Lyso PC (LPC) [M+CH3COO]-", None
            ),
        )
        self.tab_c_lipidclass_cmb.setItemText(
            16 - 1,
            QtWidgets.QApplication.translate(
                "MainWindow", "Lyso PE (LPE) [M-H]-", None
            ),
        )
        self.tab_c_lipidclass_cmb.setItemText(
            17 - 1,
            QtWidgets.QApplication.translate(
                "MainWindow", "Lyso PG (LPG) [M-H]-", None
            ),
        )
        self.tab_c_lipidclass_cmb.setItemText(
            18 - 1,
            QtWidgets.QApplication.translate(
                "MainWindow", "Lyso PI (LPI) [M-H]-", None
            ),
        )
        self.tab_c_lipidclass_cmb.setItemText(
            19 - 1,
            QtWidgets.QApplication.translate(
                "MainWindow", "Lyso PS (LPS) [M-H]-", None
            ),
        )
        self.label_17.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "MS/MS level query window:", None
            )
        )
        self.tab_c_lmms2ppm_spb.setSuffix(
            QtWidgets.QApplication.translate("MainWindow", " ppm", None)
        )
        self.tab_c_lmms2ppm_spb.setPrefix(
            QtWidgets.QApplication.translate("MainWindow", "+/- ", None)
        )
        self.tab_c_lmrun_pb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Generate Lipid Master table >>>", None
            )
        )
        self.tabframe.setTabText(
            self.tabframe.indexOf(self.cfg_tab),
            QtWidgets.QApplication.translate("MainWindow", "Settings", None),
        )
        self.label_2.setText(
            QtWidgets.QApplication.translate("MainWindow", "Lipid Class", None)
        )
        self.lipidgen_tab1_lipidclass_cmb.setItemText(
            0, QtWidgets.QApplication.translate("MainWindow", "Phospholipid", None)
        )
        self.lipidgen_tab1_lipidclass_cmb.setItemText(
            1, QtWidgets.QApplication.translate("MainWindow", "Triacylglycerol", None)
        )
        self.lipidgen_tab1_lipidclass_cmb.setItemText(
            2, QtWidgets.QApplication.translate("MainWindow", "Sphingolipid", None)
        )
        self.lipidgen_tab1_lipidclass_cmb.setItemText(
            3, QtWidgets.QApplication.translate("MainWindow", "Ceramide", None)
        )
        self.lipidgen_tab1_lipidclass_cmb.setItemText(
            4, QtWidgets.QApplication.translate("MainWindow", "Cardiolipin", None)
        )
        self.lipidgen_tab1_lipidclass_cmb.setItemText(
            5, QtWidgets.QApplication.translate("MainWindow", "Diacylglycerol", None)
        )
        self.lipidgen_tab1_all_sn_chb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Apply settings to all sn positions", None
            )
        )
        self.label_59.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Save FA whitelist as .xlsx:", None
            )
        )
        self.lipidgen_tab1_savelist_pb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Save as", None)
        )
        self.lipidgen_tab1_sn1_lb.setText(
            QtWidgets.QApplication.translate("MainWindow", "sn1", None)
        )
        self.lipidgen_tab1_sn1allfa_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Select all", None)
        )
        self.sn1_fa12x0_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA12:0", None)
        )
        self.sn1_fa22x5_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA22:5", None)
        )
        self.sn1_fa18x2_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA18:2", None)
        )
        self.sn1_fa18x1_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA18:1", None)
        )
        self.sn1_fa22x4_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA22:4", None)
        )
        self.sn1_fa20x5_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA20:5", None)
        )
        self.sn1_fa18x0_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA18:0", None)
        )
        self.sn1_fa22x6_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA22:6", None)
        )
        self.sn1_fa16x1_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA16:1", None)
        )
        self.sn1_fa14x0_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA14:0", None)
        )
        self.sn1_fa20x3_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA20:3", None)
        )
        self.sn1_fa16x0_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA16:0", None)
        )
        self.sn1_fa18x3_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA18:3", None)
        )
        self.sn1_fa20x4_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA20:4", None)
        )
        self.sn1_showop_chb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Consider O- & P- (PC,PE only)", None
            )
        )
        self.sn1_o16x0_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "O-16:0", None)
        )
        self.sn1_o18x0_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "O-18:0", None)
        )
        self.sn1_o20x0_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "O-20:0", None)
        )
        self.sn1_p16x0_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "P-16:0", None)
        )
        self.sn1_p18x0_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "P-18:0", None)
        )
        self.sn1_p20x0_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "P-20:0", None)
        )
        self.label_22.setText(
            QtWidgets.QApplication.translate("MainWindow", "Other residues*:", None)
        )
        self.label_62.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                '*: Add additional FA residue(s) using corresponging abbreviation(s) and put ", " in between e.g. FA17:0, FA14:1',
                None,
            )
        )
        self.lipidgen_tab1_sn2_lb.setText(
            QtWidgets.QApplication.translate("MainWindow", "sn2", None)
        )
        self.lipidgen_tab1_sn2allfa_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Select all", None)
        )
        self.sn2_fa22x5_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA22:5", None)
        )
        self.sn2_fa22x4_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA22:4", None)
        )
        self.sn2_fa20x5_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA20:5", None)
        )
        self.sn2_fa22x6_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA22:6", None)
        )
        self.sn2_fa18x1_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA18:1", None)
        )
        self.sn2_fa14x0_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA14:0", None)
        )
        self.sn2_fa16x0_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA16:0", None)
        )
        self.sn2_fa18x3_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA18:3", None)
        )
        self.sn2_fa18x2_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA18:2", None)
        )
        self.sn2_fa18x0_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA18:0", None)
        )
        self.sn2_fa16x1_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA16:1", None)
        )
        self.sn2_fa20x4_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA20:4", None)
        )
        self.sn2_fa20x3_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA20:3", None)
        )
        self.sn2_fa12x0_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA12:0", None)
        )
        self.label_41.setText(
            QtWidgets.QApplication.translate("MainWindow", "Other residues*:", None)
        )
        self.lipidgen_tab1_sn3_lb.setText(
            QtWidgets.QApplication.translate("MainWindow", "sn3", None)
        )
        self.lipidgen_tab1_sn3allfa_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Select all", None)
        )
        self.sn3_fa14x0_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA14:0", None)
        )
        self.sn3_fa22x4_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA22:4", None)
        )
        self.sn3_fa18x1_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA18:1", None)
        )
        self.sn3_fa22x5_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA22:5", None)
        )
        self.sn3_fa22x6_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA22:6", None)
        )
        self.sn3_fa18x3_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA18:3", None)
        )
        self.sn3_fa16x1_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA16:1", None)
        )
        self.sn3_fa20x4_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA20:4", None)
        )
        self.sn3_fa18x0_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA18:0", None)
        )
        self.sn3_fa20x5_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA20:5", None)
        )
        self.sn3_fa18x2_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA18:2", None)
        )
        self.sn3_fa16x0_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA16:0", None)
        )
        self.sn3_fa12x0_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA12:0", None)
        )
        self.sn3_fa20x3_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA20:3", None)
        )
        self.label_61.setText(
            QtWidgets.QApplication.translate("MainWindow", "Other residues*:", None)
        )
        self.lipidgen_tab1_sn4_lb.setText(
            QtWidgets.QApplication.translate("MainWindow", "sn4", None)
        )
        self.lipidgen_tab1_sn3allfa_chb_2.setText(
            QtWidgets.QApplication.translate("MainWindow", "Select all", None)
        )
        self.sn4_fa22x6_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA22:6", None)
        )
        self.sn4_fa22x5_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA22:5", None)
        )
        self.sn4_fa22x4_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA22:4", None)
        )
        self.sn4_fa20x5_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA20:5", None)
        )
        self.sn4_fa20x4_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA20:4", None)
        )
        self.sn4_fa20x3_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA20:3", None)
        )
        self.sn4_fa18x3_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA18:3", None)
        )
        self.sn4_fa18x2_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA18:2", None)
        )
        self.sn4_fa18x1_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA18:1", None)
        )
        self.sn4_fa18x0_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA18:0", None)
        )
        self.sn4_fa16x1_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA16:1", None)
        )
        self.sn4_fa16x0_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA16:0", None)
        )
        self.sn4_fa14x0_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA14:0", None)
        )
        self.sn4_fa12x0_chb.setText(
            QtWidgets.QApplication.translate("MainWindow", "FA12:0", None)
        )
        self.label_66.setText(
            QtWidgets.QApplication.translate("MainWindow", "Other residues*:", None)
        )
        self.lipidgen_tab1_genlist_pb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Generate FA whitelist ", None
            )
        )
        self.lipidgen_tabframe.setTabText(
            self.lipidgen_tabframe.indexOf(self.lipidgen_tab1),
            QtWidgets.QApplication.translate(
                "MainWindow", "Lipid Generator Wizard", None
            ),
        )
        self.label_60.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Load FA white list from .xlsx:", None
            )
        )
        self.lipidgen_tab2_loadlist_pb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Open", None)
        )
        self.lipidgen_tab2_loadlist_lb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Sucessfully loaded!", None)
        )
        self.lipidgen_tab2_savemasterlist_lb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Save lipid master list as .xlsx:", None
            )
        )
        self.lipidgen_tab2_savemasterlist_pb.setText(
            QtWidgets.QApplication.translate("MainWindow", "Save as", None)
        )
        self.lipidgen_tab1_genmasterlist_pb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Generate lipid master list >>> ", None
            )
        )
        self.lipidgen_tabframe.setTabText(
            self.lipidgen_tabframe.indexOf(self.lipidgen_tab2),
            QtWidgets.QApplication.translate(
                "MainWindow", "Create lipid master list from FA white list", None
            ),
        )
        self.tabframe.setTabText(
            self.tabframe.indexOf(self.generator_tab),
            QtWidgets.QApplication.translate(
                "MainWindow", "FA White List Generator", None
            ),
        )
        self.version_lb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                '<html><head/><body><p><span style=" font-weight:600;">LipidHunter 2 RC2 [Relased Date: 21, Jan, 2020]</span></p></body></html>',
                None,
            )
        )
        self.label_26.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "Developed by: SysMedOs team", None
            )
        )
        self.label_43.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                '<html><head/><body><p><span style=" font-weight:600;">LipidHunter is Dual-Licensed:</span></p></body></html>',
                None,
            )
        )
        self.label_42.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "       For academic & non-commercial use:", None
            )
        )
        self.label_44.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                "       For commercial use please contact the developers:",
                None,
            )
        )
        self.link_gplv2_lb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                '<html><head/><body><p>Under GPL v2 License: <a href=" https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html"><span style=" text-decoration: none; color:    #FF8C00;">https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html</span></a></p></body></html>',
                None,
            )
        )
        self.label_45.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                "<html><head/><body><p>oxlpp@bbz.uni-leipzig.de</p></body></html>",
                None,
            )
        )
        self.link_source_lb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                '<html><head/><body><p>Project page: <a href="https://github.com/SysMedOs/lipidhunter"><span style=" text-decoration: none; color:#ff8c00;">https://github.com/SysMedOs/lipidhunter</span></a></p></body></html>',
                None,
            )
        )
        self.link_tutorial_lb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                '<html><head/><body><p>Online tutorial: <a href="https://github.com/SysMedOs/lipidhunter/wiki"><span style=" text-decoration: none; color:#FF8C00;">https://github.com/SysMedOs/lipidhunter/wiki</span></a></p></body></html>',
                None,
            )
        )
        self.link_paper_lb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                '<html><head/><body><p>Please cite our publication: <a href="http://pubs.acs.org/doi/abs/10.1021/acs.analchem.7b01126"><span style=" text-decoration: none; color:#FF8C00;">Ni, Z., Angelidou, G., Lange, M., Hoffmann, R., &amp; Fedorova, M. (2017). LipidHunter identifies phospholipids by high-throughput processing of LC-MS and shotgun lipidomics datasets. Analytical chemistry, 89(17), 8800-8807.</span></a></p></body></html>',
                None,
            )
        )
        self.link_otherprojects_lb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                '<html><head/><body><p>Other SysMedOs software Projects: <a href="https://home.uni-leipzig.de/fedorova/software/"><span style=" text-decoration: none; color:#ff8c00;">https://home.uni-leipzig.de/fedorova/software/</span></a></p></body></html>',
                None,
            )
        )
        self.label_30.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                "We acknowlege all open source libraries used in LipidHunter:",
                None,
            )
        )
        self.label_31.setText(
            QtWidgets.QApplication.translate("MainWindow", "matplotlib", None)
        )
        self.link_matplotlib_lb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                '<html><head/><body><p><a href="http://www.matplotlib.org"><span style="  text-decoration: none; color:#FF8C00;">http://www.matplotlib.org</span></a><br/>John D. Hunter. Matplotlib: A 2D Graphics Environment,<br/>Computing in Science &amp; Engineering, 9, 90-95 (2007), DOI:10.1109/MCSE.2007.55 </p></body></html>',
                None,
            )
        )
        self.label_36.setText(
            QtWidgets.QApplication.translate("MainWindow", "numpy & scipy", None)
        )
        self.link_numpy_lb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                '<html><head/><body><p><a href="https://www.numpy.org"><span style=" color:#ff8c00;">https://www.numpy.org</span></a> &amp; <a href="https://www.scipy.org"><span style=" color:#ff8c00;">https://www.scipy.org</span></a><br/>Stéfan van der Walt, S. Chris Colbert and Gaël Varoquaux. The NumPy Array: A Structure for Efficient Numerical Computation, <br/>Computing in Science &amp; Engineering, 13, 22-30 (2011), DOI:10.1109/MCSE.2011.37 </p></body></html>',
                None,
            )
        )
        self.label_32.setText(
            QtWidgets.QApplication.translate("MainWindow", "pandas", None)
        )
        self.link_pandas_lb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                '<html><head/><body><p><a href="http://pandas.pydata.org"><span style="  text-decoration: none; color:#FF8C00;">http://pandas.pydata.org</span></a><br/>Wes McKinney. Data Structures for Statistical Computing in Python,<br/>Proceedings of the 9th Python in Science Conference, 51-56 (2010) </p></body></html>',
                None,
            )
        )
        self.label_33.setText(
            QtWidgets.QApplication.translate("MainWindow", "pymzML", None)
        )
        self.link_pymzml_lb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                '<html><head/><body><p><a href="https://github.com/pymzml/pymzML"><span style=" color:#ff8c00;">https://github.com/pymzml/pymzML</span></a><br/>Bald, T., Barth, J., Niehues, A., Specht, M., Hippler, M., and Fufezan, C. (2012) <br/>pymzML - Python module for high throughput bioinformatics on mass spectrometry data,<br/>Bioinformatics, DOI: 10.1093/bioinformatics/bts066 </p></body></html>',
                None,
            )
        )
        self.label_34.setText(
            QtWidgets.QApplication.translate("MainWindow", "PySide2", None)
        )
        self.link_pyside_lb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                '<html><head/><body><p><a href="https://wiki.qt.io/Qt_for_Python"><span style="  text-decoration: none; color:#FF8C00;">https://wiki.qt.io/Qt_for_Python</span></a></p></body></html>',
                None,
            )
        )
        self.label_37.setText(
            QtWidgets.QApplication.translate("MainWindow", "numba", None)
        )
        self.label_97.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                '<html><head/><body><p><a href="http://numba.pydata.org"><span style="  text-decoration: none; color:#FF8C00;">http://numba.pydata.org</span></a><br/>Siu Kwan Lam, Antoine Pitrou, and Stanley Seibert.  Numba: a LLVM-based Python JIT compiler,<br/> DOI:10.1145/2833157.2833162 </p></body></html>',
                None,
            )
        )
        self.label_46.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                "We acknowlege all projects that support the development of LipidHunter:",
                None,
            )
        )
        self.label_47.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                "BMBF - Federal Ministry of Education and Research Germany:",
                None,
            )
        )
        self.label_49.setText(
            QtWidgets.QApplication.translate(
                "MainWindow", "e:Med Systems Medicine Network:", None
            )
        )
        self.link_bmbf_lb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                '<html><head/><body><p><a href="https://www.bmbf.de/en/"><span style="  text-decoration: none; color:#FF8C00;">https://www.bmbf.de/en/</span></a></p></body></html>',
                None,
            )
        )
        self.link_emed_lb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                '<html><head/><body><p><a href="http://www.sys-med.de/en/"><span style="   text-decoration: none; color:#FF8C00;">http://www.sys-med.de/en/</span></a></p></body></html>',
                None,
            )
        )
        self.link_sysmedos_lb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                '<html><head/><body><p><a href="https://home.uni-leipzig.de/fedorova/sysmedos/"><span style="   text-decoration: none; color:#FF8C00;">https://home.uni-leipzig.de/fedorova/sysmedos/</span></a></p></body></html>',
                None,
            )
        )
        self.label_58.setText(
            QtWidgets.QApplication.translate("MainWindow", "SysMedOS Project : ", None)
        )
        self.tabframe.setTabText(
            self.tabframe.indexOf(self.tab),
            QtWidgets.QApplication.translate("MainWindow", "About LipidHunter 2", None),
        )
        self.link_uni_lb.setText(
            QtWidgets.QApplication.translate(
                "MainWindow",
                '<html><head/><body><p>LipidHunter (C) 2016-2020 | <a href="http://www.zv.uni-leipzig.de/"><span style=" color:#ff8c00;">University of Leipzig</span></a> | <a href="http://research.uni-leipzig.de/bioanalytik/"><span style=" color:#ff8c00;">AG Bioanalytik</span></a> | <a href="https://home.uni-leipzig.de/fedorova/"><span style=" color:#ff8c00;">Fedorova Research Group</span></a></p></body></html>',
                None,
            )
        )

    import LibLipidHunter.LipidHunter_rcc
