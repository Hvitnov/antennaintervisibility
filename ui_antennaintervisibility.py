# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui_antennaintervisibility.ui'
#
# Created: Mon Jun 20 16:20:16 2016
#      by: PyQt4 UI code generator 4.10.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_AntennaIntervisibility(object):
    def setupUi(self, AntennaIntervisibility):
        AntennaIntervisibility.setObjectName(_fromUtf8("AntennaIntervisibility"))
        AntennaIntervisibility.resize(434, 683)
        self.tabWidget = QtGui.QTabWidget(AntennaIntervisibility)
        self.tabWidget.setGeometry(QtCore.QRect(0, 10, 431, 671))
        self.tabWidget.setObjectName(_fromUtf8("tabWidget"))
        self.General_tab = QtGui.QWidget()
        self.General_tab.setObjectName(_fromUtf8("General_tab"))
        self.layoutWidget = QtGui.QWidget(self.General_tab)
        self.layoutWidget.setGeometry(QtCore.QRect(8, 10, 411, 618))
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.verticalLayout_3 = QtGui.QVBoxLayout(self.layoutWidget)
        self.verticalLayout_3.setMargin(0)
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.DEM_label_2 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.DEM_label_2.sizePolicy().hasHeightForWidth())
        self.DEM_label_2.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(11)
        font.setBold(True)
        font.setWeight(75)
        self.DEM_label_2.setFont(font)
        self.DEM_label_2.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.DEM_label_2.setObjectName(_fromUtf8("DEM_label_2"))
        self.verticalLayout_2.addWidget(self.DEM_label_2)
        self.horizontalLayout_5 = QtGui.QHBoxLayout()
        self.horizontalLayout_5.setObjectName(_fromUtf8("horizontalLayout_5"))
        self.label_6 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_6.sizePolicy().hasHeightForWidth())
        self.label_6.setSizePolicy(sizePolicy)
        self.label_6.setMinimumSize(QtCore.QSize(100, 0))
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.horizontalLayout_5.addWidget(self.label_6)
        self.horizontalLayout_6 = QtGui.QHBoxLayout()
        self.horizontalLayout_6.setObjectName(_fromUtf8("horizontalLayout_6"))
        self.cmbPoints = QtGui.QComboBox(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.cmbPoints.sizePolicy().hasHeightForWidth())
        self.cmbPoints.setSizePolicy(sizePolicy)
        self.cmbPoints.setEditable(False)
        self.cmbPoints.setObjectName(_fromUtf8("cmbPoints"))
        self.horizontalLayout_6.addWidget(self.cmbPoints)
        self.horizontalLayout_5.addLayout(self.horizontalLayout_6)
        self.verticalLayout_2.addLayout(self.horizontalLayout_5)
        self.verticalLayout_3.addLayout(self.verticalLayout_2)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.label_3 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_3.sizePolicy().hasHeightForWidth())
        self.label_3.setSizePolicy(sizePolicy)
        self.label_3.setMinimumSize(QtCore.QSize(100, 0))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.horizontalLayout.addWidget(self.label_3)
        self.txtObserver = QtGui.QLineEdit(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.txtObserver.sizePolicy().hasHeightForWidth())
        self.txtObserver.setSizePolicy(sizePolicy)
        self.txtObserver.setObjectName(_fromUtf8("txtObserver"))
        self.horizontalLayout.addWidget(self.txtObserver)
        self.label_10 = QtGui.QLabel(self.layoutWidget)
        self.label_10.setObjectName(_fromUtf8("label_10"))
        self.horizontalLayout.addWidget(self.label_10)
        self.cmbObsField = QtGui.QComboBox(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.cmbObsField.sizePolicy().hasHeightForWidth())
        self.cmbObsField.setSizePolicy(sizePolicy)
        self.cmbObsField.setEditable(False)
        self.cmbObsField.setObjectName(_fromUtf8("cmbObsField"))
        self.horizontalLayout.addWidget(self.cmbObsField)
        self.verticalLayout_3.addLayout(self.horizontalLayout)
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.verticalLayout_5 = QtGui.QVBoxLayout()
        self.verticalLayout_5.setObjectName(_fromUtf8("verticalLayout_5"))
        self.DEM_label_3 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.DEM_label_3.sizePolicy().hasHeightForWidth())
        self.DEM_label_3.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(11)
        font.setBold(True)
        font.setWeight(75)
        self.DEM_label_3.setFont(font)
        self.DEM_label_3.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.DEM_label_3.setObjectName(_fromUtf8("DEM_label_3"))
        self.verticalLayout_5.addWidget(self.DEM_label_3)
        self.horizontalLayout_7 = QtGui.QHBoxLayout()
        self.horizontalLayout_7.setObjectName(_fromUtf8("horizontalLayout_7"))
        self.label_9 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_9.sizePolicy().hasHeightForWidth())
        self.label_9.setSizePolicy(sizePolicy)
        self.label_9.setMinimumSize(QtCore.QSize(100, 0))
        self.label_9.setObjectName(_fromUtf8("label_9"))
        self.horizontalLayout_7.addWidget(self.label_9)
        self.horizontalLayout_8 = QtGui.QHBoxLayout()
        self.horizontalLayout_8.setObjectName(_fromUtf8("horizontalLayout_8"))
        self.cmbPointsTarget = QtGui.QComboBox(self.layoutWidget)
        self.cmbPointsTarget.setEnabled(True)
        self.cmbPointsTarget.setEditable(False)
        self.cmbPointsTarget.setObjectName(_fromUtf8("cmbPointsTarget"))
        self.horizontalLayout_8.addWidget(self.cmbPointsTarget)
        self.horizontalLayout_7.addLayout(self.horizontalLayout_8)
        self.verticalLayout_5.addLayout(self.horizontalLayout_7)
        self.horizontalLayout_9 = QtGui.QHBoxLayout()
        self.horizontalLayout_9.setObjectName(_fromUtf8("horizontalLayout_9"))
        self.label_4 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_4.sizePolicy().hasHeightForWidth())
        self.label_4.setSizePolicy(sizePolicy)
        self.label_4.setMinimumSize(QtCore.QSize(100, 0))
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.horizontalLayout_9.addWidget(self.label_4)
        self.txtTarget = QtGui.QLineEdit(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.txtTarget.sizePolicy().hasHeightForWidth())
        self.txtTarget.setSizePolicy(sizePolicy)
        self.txtTarget.setObjectName(_fromUtf8("txtTarget"))
        self.horizontalLayout_9.addWidget(self.txtTarget)
        self.label_11 = QtGui.QLabel(self.layoutWidget)
        self.label_11.setObjectName(_fromUtf8("label_11"))
        self.horizontalLayout_9.addWidget(self.label_11)
        self.cmbTargetField = QtGui.QComboBox(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.cmbTargetField.sizePolicy().hasHeightForWidth())
        self.cmbTargetField.setSizePolicy(sizePolicy)
        self.cmbTargetField.setEditable(False)
        self.cmbTargetField.setObjectName(_fromUtf8("cmbTargetField"))
        self.horizontalLayout_9.addWidget(self.cmbTargetField)
        self.verticalLayout_5.addLayout(self.horizontalLayout_9)
        self.verticalLayout.addLayout(self.verticalLayout_5)
        self.DEM_label = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.DEM_label.sizePolicy().hasHeightForWidth())
        self.DEM_label.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(11)
        font.setBold(True)
        font.setWeight(75)
        self.DEM_label.setFont(font)
        self.DEM_label.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.DEM_label.setObjectName(_fromUtf8("DEM_label"))
        self.verticalLayout.addWidget(self.DEM_label)
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.label = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        self.label.setMinimumSize(QtCore.QSize(100, 0))
        self.label.setObjectName(_fromUtf8("label"))
        self.horizontalLayout_2.addWidget(self.label)
        self.horizontalLayout_4 = QtGui.QHBoxLayout()
        self.horizontalLayout_4.setObjectName(_fromUtf8("horizontalLayout_4"))
        self.cmbRaster = QtGui.QComboBox(self.layoutWidget)
        self.cmbRaster.setObjectName(_fromUtf8("cmbRaster"))
        self.horizontalLayout_4.addWidget(self.cmbRaster)
        self.horizontalLayout_2.addLayout(self.horizontalLayout_4)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.verticalLayout_6 = QtGui.QVBoxLayout()
        self.verticalLayout_6.setObjectName(_fromUtf8("verticalLayout_6"))
        self.DEM_label_5 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.DEM_label_5.sizePolicy().hasHeightForWidth())
        self.DEM_label_5.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(11)
        font.setBold(True)
        font.setWeight(75)
        self.DEM_label_5.setFont(font)
        self.DEM_label_5.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.DEM_label_5.setObjectName(_fromUtf8("DEM_label_5"))
        self.verticalLayout_6.addWidget(self.DEM_label_5)
        self.horizontalLayout_19 = QtGui.QHBoxLayout()
        self.horizontalLayout_19.setSpacing(6)
        self.horizontalLayout_19.setObjectName(_fromUtf8("horizontalLayout_19"))
        self.horizontalLayout_20 = QtGui.QHBoxLayout()
        self.horizontalLayout_20.setObjectName(_fromUtf8("horizontalLayout_20"))
        self.cmdBrowse = QtGui.QPushButton(self.layoutWidget)
        self.cmdBrowse.setObjectName(_fromUtf8("cmdBrowse"))
        self.horizontalLayout_20.addWidget(self.cmdBrowse)
        self.txtOutput = QtGui.QPlainTextEdit(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Ignored)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.txtOutput.sizePolicy().hasHeightForWidth())
        self.txtOutput.setSizePolicy(sizePolicy)
        self.txtOutput.setMaximumSize(QtCore.QSize(16777215, 30))
        self.txtOutput.setObjectName(_fromUtf8("txtOutput"))
        self.horizontalLayout_20.addWidget(self.txtOutput)
        self.horizontalLayout_19.addLayout(self.horizontalLayout_20)
        self.verticalLayout_6.addLayout(self.horizontalLayout_19)
        self.DEM_label_7 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.DEM_label_7.sizePolicy().hasHeightForWidth())
        self.DEM_label_7.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(11)
        font.setBold(True)
        font.setWeight(75)
        self.DEM_label_7.setFont(font)
        self.DEM_label_7.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.DEM_label_7.setObjectName(_fromUtf8("DEM_label_7"))
        self.verticalLayout_6.addWidget(self.DEM_label_7)
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        self.label_8 = QtGui.QLabel(self.layoutWidget)
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.horizontalLayout_3.addWidget(self.label_8)
        self.txtRadius = QtGui.QLineEdit(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.txtRadius.sizePolicy().hasHeightForWidth())
        self.txtRadius.setSizePolicy(sizePolicy)
        self.txtRadius.setObjectName(_fromUtf8("txtRadius"))
        self.horizontalLayout_3.addWidget(self.txtRadius)
        self.verticalLayout_6.addLayout(self.horizontalLayout_3)
        self.horizontalLayout_10 = QtGui.QHBoxLayout()
        self.horizontalLayout_10.setObjectName(_fromUtf8("horizontalLayout_10"))
        self.checkBox_Fresnel = QtGui.QRadioButton(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.checkBox_Fresnel.sizePolicy().hasHeightForWidth())
        self.checkBox_Fresnel.setSizePolicy(sizePolicy)
        self.checkBox_Fresnel.setChecked(False)
        self.checkBox_Fresnel.setObjectName(_fromUtf8("checkBox_Fresnel"))
        self.horizontalLayout_10.addWidget(self.checkBox_Fresnel)
        self.verticalLayout_6.addLayout(self.horizontalLayout_10)
        self.verticalLayout_4 = QtGui.QVBoxLayout()
        self.verticalLayout_4.setObjectName(_fromUtf8("verticalLayout_4"))
        self.DEM_label_6 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.DEM_label_6.sizePolicy().hasHeightForWidth())
        self.DEM_label_6.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(11)
        font.setBold(True)
        font.setWeight(75)
        self.DEM_label_6.setFont(font)
        self.DEM_label_6.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.DEM_label_6.setObjectName(_fromUtf8("DEM_label_6"))
        self.verticalLayout_4.addWidget(self.DEM_label_6)
        self.chkIntervisibility = QtGui.QRadioButton(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.chkIntervisibility.sizePolicy().hasHeightForWidth())
        self.chkIntervisibility.setSizePolicy(sizePolicy)
        self.chkIntervisibility.setChecked(True)
        self.chkIntervisibility.setObjectName(_fromUtf8("chkIntervisibility"))
        self.verticalLayout_4.addWidget(self.chkIntervisibility)
        self.chkIntervisibility_2 = QtGui.QRadioButton(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.chkIntervisibility_2.sizePolicy().hasHeightForWidth())
        self.chkIntervisibility_2.setSizePolicy(sizePolicy)
        self.chkIntervisibility_2.setChecked(False)
        self.chkIntervisibility_2.setObjectName(_fromUtf8("chkIntervisibility_2"))
        self.verticalLayout_4.addWidget(self.chkIntervisibility_2)
        self.verticalLayout_6.addLayout(self.verticalLayout_4)
        self.horizontalLayout_14 = QtGui.QHBoxLayout()
        self.horizontalLayout_14.setObjectName(_fromUtf8("horizontalLayout_14"))
        self.label_13 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_13.sizePolicy().hasHeightForWidth())
        self.label_13.setSizePolicy(sizePolicy)
        self.label_13.setMinimumSize(QtCore.QSize(100, 0))
        self.label_13.setObjectName(_fromUtf8("label_13"))
        self.horizontalLayout_14.addWidget(self.label_13)
        self.OutputTypeSelector = QtGui.QComboBox(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.OutputTypeSelector.sizePolicy().hasHeightForWidth())
        self.OutputTypeSelector.setSizePolicy(sizePolicy)
        self.OutputTypeSelector.setEditable(False)
        self.OutputTypeSelector.setObjectName(_fromUtf8("OutputTypeSelector"))
        self.OutputTypeSelector.addItem(_fromUtf8(""))
        self.OutputTypeSelector.addItem(_fromUtf8(""))
        self.horizontalLayout_14.addWidget(self.OutputTypeSelector)
        self.verticalLayout_6.addLayout(self.horizontalLayout_14)
        self.verticalLayout.addLayout(self.verticalLayout_6)
        self.verticalLayout_3.addLayout(self.verticalLayout)
        self.horizontalLayout_13 = QtGui.QHBoxLayout()
        self.horizontalLayout_13.setObjectName(_fromUtf8("horizontalLayout_13"))
        self.buttonBox = QtGui.QDialogButtonBox(self.layoutWidget)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.horizontalLayout_13.addWidget(self.buttonBox)
        self.verticalLayout_3.addLayout(self.horizontalLayout_13)
        self.tabWidget.addTab(self.General_tab, _fromUtf8(""))
        self.tab = QtGui.QWidget()
        self.tab.setObjectName(_fromUtf8("tab"))
        self.textBrowser_2 = QtGui.QTextBrowser(self.tab)
        self.textBrowser_2.setGeometry(QtCore.QRect(0, 0, 431, 641))
        self.textBrowser_2.setObjectName(_fromUtf8("textBrowser_2"))
        self.tabWidget.addTab(self.tab, _fromUtf8(""))
        self.About_tab = QtGui.QWidget()
        self.About_tab.setObjectName(_fromUtf8("About_tab"))
        self.textBrowser = QtGui.QTextBrowser(self.About_tab)
        self.textBrowser.setGeometry(QtCore.QRect(0, 0, 431, 641))
        self.textBrowser.setFocusPolicy(QtCore.Qt.WheelFocus)
        self.textBrowser.setObjectName(_fromUtf8("textBrowser"))
        self.tabWidget.addTab(self.About_tab, _fromUtf8(""))

        self.retranslateUi(AntennaIntervisibility)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), AntennaIntervisibility.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), AntennaIntervisibility.reject)
        QtCore.QMetaObject.connectSlotsByName(AntennaIntervisibility)

    def retranslateUi(self, AntennaIntervisibility):
        AntennaIntervisibility.setWindowTitle(_translate("AntennaIntervisibility", "Antenna Intervisibility", None))
        self.DEM_label_2.setText(_translate("AntennaIntervisibility", "Observation Points", None))
        self.label_6.setText(_translate("AntennaIntervisibility", "Select Map Layer", None))
        self.label_3.setText(_translate("AntennaIntervisibility", "Observer height", None))
        self.txtObserver.setText(_translate("AntennaIntervisibility", "1.6", None))
        self.label_10.setText(_translate("AntennaIntervisibility", "or field:", None))
        self.DEM_label_3.setText(_translate("AntennaIntervisibility", "Target Points", None))
        self.label_9.setText(_translate("AntennaIntervisibility", "Select Map Layer", None))
        self.label_4.setText(_translate("AntennaIntervisibility", "Target height", None))
        self.txtTarget.setText(_translate("AntennaIntervisibility", "0", None))
        self.label_11.setText(_translate("AntennaIntervisibility", "or field:", None))
        self.DEM_label.setText(_translate("AntennaIntervisibility", "Digital Elevation Model (Raster)", None))
        self.label.setText(_translate("AntennaIntervisibility", "Select Map Layer", None))
        self.DEM_label_5.setText(_translate("AntennaIntervisibility", "Output File", None))
        self.cmdBrowse.setText(_translate("AntennaIntervisibility", "Browse", None))
        self.txtOutput.setPlainText(_translate("AntennaIntervisibility", "/home/hvitnov/test", None))
        self.DEM_label_7.setText(_translate("AntennaIntervisibility", "Settings", None))
        self.label_8.setText(_translate("AntennaIntervisibility", "Search radius", None))
        self.txtRadius.setText(_translate("AntennaIntervisibility", "2000", None))
        self.checkBox_Fresnel.setText(_translate("AntennaIntervisibility", "Enable Fresnel Zones", None))
        self.DEM_label_6.setText(_translate("AntennaIntervisibility", "Algorithm Type", None))
        self.chkIntervisibility.setText(_translate("AntennaIntervisibility", "Intervisibility - Zoran c version", None))
        self.chkIntervisibility_2.setText(_translate("AntennaIntervisibility", "Intervisibility - Something else", None))
        self.label_13.setText(_translate("AntennaIntervisibility", "Output Type:", None))
        self.OutputTypeSelector.setItemText(0, _translate("AntennaIntervisibility", "Highest Points in Polygon", None))
        self.OutputTypeSelector.setItemText(1, _translate("AntennaIntervisibility", "Lowest Points in Polygon", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.General_tab), _translate("AntennaIntervisibility", "Setup", None))
        self.textBrowser_2.setHtml(_translate("AntennaIntervisibility", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Noto Sans\'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Ubuntu\'; font-size:10pt; font-weight:600;\">Observer Points:</span><span style=\" font-family:\'Ubuntu\'; font-size:10pt;\"> A Vector layer containing points from which to measure intervisibility to targets</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:\'Ubuntu\'; font-size:10pt;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Ubuntu\'; font-size:10pt; font-weight:600;\">Target Points:</span><span style=\" font-family:\'Ubuntu\'; font-size:10pt;\"> A Vector layer containing points from which to measure intervisibility to observers</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:\'Ubuntu\'; font-size:10pt;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Ubuntu\'; font-size:10pt; font-weight:600;\">Digital Elevation Model:</span><span style=\" font-family:\'Ubuntu\'; font-size:10pt;\"> A raster layer representing a digital elevation model</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:\'Ubuntu\'; font-size:10pt;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Ubuntu\'; font-size:10pt; font-weight:600;\">Options:</span><span style=\" font-family:\'Ubuntu\'; font-size:10pt;\"> </span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Ubuntu\'; font-size:10pt;\">Radius: The furthest distance visible from observer (in map units)</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Ubuntu\'; font-size:10pt;\">Enable Fresnel Zones: Consider Fresnel (Radio propagation) Zones when meausuring intervisibility</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Ubuntu\'; font-size:10pt; font-weight:600;\">Output file:</span><span style=\" font-family:\'Ubuntu\'; font-size:10pt;\">  A .shp file containing a field for highest as well as lowest points in polygons</span></p></body></html>", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("AntennaIntervisibility", "Help", None))
        self.textBrowser.setHtml(_translate("AntennaIntervisibility", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Noto Sans\'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:\'Ubuntu\'; font-size:11pt;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Ubuntu\'; font-size:12pt; font-weight:600;\">Version 1.0</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Ubuntu\'; font-size:11pt;\">Developed by Jakob Hvitnov, Roskilde Universitetscenter, Denmark</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:\'Ubuntu\'; font-size:11pt;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Ubuntu\'; font-size:11pt;\">Intervisibility Algorithm 1 from &quot;ViewshedAnalysis&quot; plugin by Zoran Čučković</span></p></body></html>", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.About_tab), _translate("AntennaIntervisibility", "About", None))

