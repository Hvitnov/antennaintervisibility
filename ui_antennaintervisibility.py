# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui_antennaintervisibility.ui'
#
# Created: Sat Apr 23 16:44:57 2016
#      by: PyQt4 UI code generator 4.10.2
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
        AntennaIntervisibility.resize(448, 622)
        self.tabWidget = QtGui.QTabWidget(AntennaIntervisibility)
        self.tabWidget.setGeometry(QtCore.QRect(10, 10, 427, 571))
        self.tabWidget.setObjectName(_fromUtf8("tabWidget"))
        self.General_tab = QtGui.QWidget()
        self.General_tab.setObjectName(_fromUtf8("General_tab"))
        self.cmbPoints = QtGui.QComboBox(self.General_tab)
        self.cmbPoints.setGeometry(QtCore.QRect(10, 90, 201, 22))
        self.cmbPoints.setEditable(False)
        self.cmbPoints.setObjectName(_fromUtf8("cmbPoints"))
        self.label = QtGui.QLabel(self.General_tab)
        self.label.setGeometry(QtCore.QRect(10, 20, 111, 16))
        self.label.setObjectName(_fromUtf8("label"))
        self.label_2 = QtGui.QLabel(self.General_tab)
        self.label_2.setGeometry(QtCore.QRect(10, 70, 131, 16))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.cmbRaster = QtGui.QComboBox(self.General_tab)
        self.cmbRaster.setGeometry(QtCore.QRect(10, 40, 201, 22))
        self.cmbRaster.setObjectName(_fromUtf8("cmbRaster"))
        self.txtOutput = QtGui.QPlainTextEdit(self.General_tab)
        self.txtOutput.setGeometry(QtCore.QRect(220, 40, 181, 41))
        self.txtOutput.setObjectName(_fromUtf8("txtOutput"))
        self.label_5 = QtGui.QLabel(self.General_tab)
        self.label_5.setGeometry(QtCore.QRect(220, 20, 91, 16))
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.cmdBrowse = QtGui.QPushButton(self.General_tab)
        self.cmdBrowse.setGeometry(QtCore.QRect(320, 90, 75, 31))
        self.cmdBrowse.setObjectName(_fromUtf8("cmdBrowse"))
        self.Output_groupBox = QtGui.QGroupBox(self.General_tab)
        self.Output_groupBox.setGeometry(QtCore.QRect(10, 350, 391, 91))
        self.Output_groupBox.setObjectName(_fromUtf8("Output_groupBox"))
        self.chkBinary = QtGui.QRadioButton(self.Output_groupBox)
        self.chkBinary.setGeometry(QtCore.QRect(10, 30, 141, 22))
        self.chkBinary.setObjectName(_fromUtf8("chkBinary"))
        self.chkInvisibility = QtGui.QRadioButton(self.Output_groupBox)
        self.chkInvisibility.setGeometry(QtCore.QRect(200, 30, 141, 22))
        self.chkInvisibility.setObjectName(_fromUtf8("chkInvisibility"))
        self.chkIntervisibility = QtGui.QRadioButton(self.Output_groupBox)
        self.chkIntervisibility.setGeometry(QtCore.QRect(10, 60, 161, 17))
        self.chkIntervisibility.setObjectName(_fromUtf8("chkIntervisibility"))
        self.chkHorizon = QtGui.QRadioButton(self.Output_groupBox)
        self.chkHorizon.setGeometry(QtCore.QRect(200, 60, 181, 16))
        self.chkHorizon.setObjectName(_fromUtf8("chkHorizon"))
        self.groupBox_3 = QtGui.QGroupBox(self.General_tab)
        self.groupBox_3.setGeometry(QtCore.QRect(10, 170, 391, 167))
        self.groupBox_3.setObjectName(_fromUtf8("groupBox_3"))
        self.label_16 = QtGui.QLabel(self.groupBox_3)
        self.label_16.setGeometry(QtCore.QRect(260, 140, 121, 20))
        self.label_16.setObjectName(_fromUtf8("label_16"))
        self.txtAdaptRadiusObs = QtGui.QLineEdit(self.groupBox_3)
        self.txtAdaptRadiusObs.setGeometry(QtCore.QRect(40, 140, 31, 20))
        self.txtAdaptRadiusObs.setObjectName(_fromUtf8("txtAdaptRadiusObs"))
        self.label_15 = QtGui.QLabel(self.groupBox_3)
        self.label_15.setGeometry(QtCore.QRect(80, 140, 131, 20))
        self.label_15.setObjectName(_fromUtf8("label_15"))
        self.txtAdaptRadiusTarget = QtGui.QLineEdit(self.groupBox_3)
        self.txtAdaptRadiusTarget.setGeometry(QtCore.QRect(220, 140, 31, 20))
        self.txtAdaptRadiusTarget.setObjectName(_fromUtf8("txtAdaptRadiusTarget"))
        self.cmbTargetField = QtGui.QComboBox(self.groupBox_3)
        self.cmbTargetField.setGeometry(QtCore.QRect(260, 90, 121, 22))
        self.cmbTargetField.setEditable(False)
        self.cmbTargetField.setObjectName(_fromUtf8("cmbTargetField"))
        self.label_10 = QtGui.QLabel(self.groupBox_3)
        self.label_10.setGeometry(QtCore.QRect(190, 60, 61, 20))
        self.label_10.setObjectName(_fromUtf8("label_10"))
        self.txtObserver = QtGui.QLineEdit(self.groupBox_3)
        self.txtObserver.setGeometry(QtCore.QRect(130, 60, 51, 20))
        self.txtObserver.setObjectName(_fromUtf8("txtObserver"))
        self.cmbObsField = QtGui.QComboBox(self.groupBox_3)
        self.cmbObsField.setGeometry(QtCore.QRect(260, 60, 121, 22))
        self.cmbObsField.setEditable(False)
        self.cmbObsField.setObjectName(_fromUtf8("cmbObsField"))
        self.label_11 = QtGui.QLabel(self.groupBox_3)
        self.label_11.setGeometry(QtCore.QRect(190, 90, 61, 20))
        self.label_11.setObjectName(_fromUtf8("label_11"))
        self.label_4 = QtGui.QLabel(self.groupBox_3)
        self.label_4.setGeometry(QtCore.QRect(10, 90, 101, 16))
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.label_3 = QtGui.QLabel(self.groupBox_3)
        self.label_3.setGeometry(QtCore.QRect(10, 60, 121, 16))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.txtTarget = QtGui.QLineEdit(self.groupBox_3)
        self.txtTarget.setGeometry(QtCore.QRect(130, 90, 51, 20))
        self.txtTarget.setObjectName(_fromUtf8("txtTarget"))
        self.label_6 = QtGui.QLabel(self.groupBox_3)
        self.label_6.setGeometry(QtCore.QRect(10, 120, 271, 16))
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.txtRadius = QtGui.QLineEdit(self.groupBox_3)
        self.txtRadius.setGeometry(QtCore.QRect(130, 30, 51, 21))
        self.txtRadius.setObjectName(_fromUtf8("txtRadius"))
        self.label_8 = QtGui.QLabel(self.groupBox_3)
        self.label_8.setGeometry(QtCore.QRect(10, 30, 91, 16))
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.groupBox = QtGui.QGroupBox(self.General_tab)
        self.groupBox.setGeometry(QtCore.QRect(10, 450, 391, 81))
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.chkCurvature = QtGui.QCheckBox(self.groupBox)
        self.chkCurvature.setGeometry(QtCore.QRect(20, 50, 161, 22))
        self.chkCurvature.setObjectName(_fromUtf8("chkCurvature"))
        self.txtRefraction = QtGui.QLineEdit(self.groupBox)
        self.txtRefraction.setEnabled(False)
        self.txtRefraction.setGeometry(QtCore.QRect(190, 50, 51, 21))
        self.txtRefraction.setObjectName(_fromUtf8("txtRefraction"))
        self.label_17 = QtGui.QLabel(self.groupBox)
        self.label_17.setGeometry(QtCore.QRect(250, 50, 161, 20))
        self.label_17.setObjectName(_fromUtf8("label_17"))
        self.chkCumulative = QtGui.QCheckBox(self.groupBox)
        self.chkCumulative.setEnabled(True)
        self.chkCumulative.setGeometry(QtCore.QRect(20, 20, 261, 21))
        self.chkCumulative.setObjectName(_fromUtf8("chkCumulative"))
        self.cmbPointsTarget = QtGui.QComboBox(self.General_tab)
        self.cmbPointsTarget.setEnabled(True)
        self.cmbPointsTarget.setGeometry(QtCore.QRect(10, 140, 201, 22))
        self.cmbPointsTarget.setEditable(False)
        self.cmbPointsTarget.setObjectName(_fromUtf8("cmbPointsTarget"))
        self.label_7 = QtGui.QLabel(self.General_tab)
        self.label_7.setGeometry(QtCore.QRect(10, 120, 191, 20))
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.tabWidget.addTab(self.General_tab, _fromUtf8(""))
        self.tab = QtGui.QWidget()
        self.tab.setObjectName(_fromUtf8("tab"))
        self.textBrowser_2 = QtGui.QTextBrowser(self.tab)
        self.textBrowser_2.setGeometry(QtCore.QRect(0, 0, 421, 501))
        self.textBrowser_2.setObjectName(_fromUtf8("textBrowser_2"))
        self.label_19 = QtGui.QLabel(self.tab)
        self.label_19.setGeometry(QtCore.QRect(50, 510, 301, 21))
        self.label_19.setTextFormat(QtCore.Qt.RichText)
        self.label_19.setOpenExternalLinks(True)
        self.label_19.setObjectName(_fromUtf8("label_19"))
        self.tabWidget.addTab(self.tab, _fromUtf8(""))
        self.About_tab = QtGui.QWidget()
        self.About_tab.setObjectName(_fromUtf8("About_tab"))
        self.textBrowser = QtGui.QTextBrowser(self.About_tab)
        self.textBrowser.setGeometry(QtCore.QRect(0, 0, 411, 501))
        self.textBrowser.setFocusPolicy(QtCore.Qt.WheelFocus)
        self.textBrowser.setObjectName(_fromUtf8("textBrowser"))
        self.label_20 = QtGui.QLabel(self.About_tab)
        self.label_20.setGeometry(QtCore.QRect(50, 510, 301, 21))
        self.label_20.setTextFormat(QtCore.Qt.RichText)
        self.label_20.setOpenExternalLinks(True)
        self.label_20.setObjectName(_fromUtf8("label_20"))
        self.tabWidget.addTab(self.About_tab, _fromUtf8(""))
        self.buttonBox = QtGui.QDialogButtonBox(AntennaIntervisibility)
        self.buttonBox.setGeometry(QtCore.QRect(270, 590, 161, 32))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))

        self.retranslateUi(AntennaIntervisibility)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), AntennaIntervisibility.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), AntennaIntervisibility.reject)
        QtCore.QObject.connect(self.chkCurvature, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.txtRefraction.setEnabled)
        QtCore.QMetaObject.connectSlotsByName(AntennaIntervisibility)

    def retranslateUi(self, AntennaIntervisibility):
        AntennaIntervisibility.setWindowTitle(_translate("AntennaIntervisibility", "Advanced viewshed analysis", None))
        self.label.setText(_translate("AntennaIntervisibility", "Elevation raster", None))
        self.label_2.setText(_translate("AntennaIntervisibility", "Observation points", None))
        self.label_5.setText(_translate("AntennaIntervisibility", "Output file", None))
        self.cmdBrowse.setText(_translate("AntennaIntervisibility", "Browse", None))
        self.Output_groupBox.setTitle(_translate("AntennaIntervisibility", "Output ", None))
        self.chkBinary.setText(_translate("AntennaIntervisibility", "Binary viewshed", None))
        self.chkInvisibility.setText(_translate("AntennaIntervisibility", "Invisibility depth", None))
        self.chkIntervisibility.setText(_translate("AntennaIntervisibility", "Intervisibility", None))
        self.chkHorizon.setText(_translate("AntennaIntervisibility", "Horizon", None))
        self.groupBox_3.setTitle(_translate("AntennaIntervisibility", "Settings ", None))
        self.label_16.setText(_translate("AntennaIntervisibility", "pixels for target", None))
        self.txtAdaptRadiusObs.setText(_translate("AntennaIntervisibility", "0", None))
        self.label_15.setText(_translate("AntennaIntervisibility", "pixels for observer", None))
        self.txtAdaptRadiusTarget.setText(_translate("AntennaIntervisibility", "0", None))
        self.label_10.setText(_translate("AntennaIntervisibility", "or field:", None))
        self.txtObserver.setText(_translate("AntennaIntervisibility", "1.6", None))
        self.label_11.setText(_translate("AntennaIntervisibility", "or field:", None))
        self.label_4.setText(_translate("AntennaIntervisibility", "Target height", None))
        self.label_3.setText(_translate("AntennaIntervisibility", "Observer height", None))
        self.txtTarget.setText(_translate("AntennaIntervisibility", "0", None))
        self.label_6.setText(_translate("AntennaIntervisibility", "Adapt to highest point at a distance of: ", None))
        self.txtRadius.setText(_translate("AntennaIntervisibility", "5000", None))
        self.label_8.setText(_translate("AntennaIntervisibility", "Search radius", None))
        self.groupBox.setTitle(_translate("AntennaIntervisibility", "Options", None))
        self.chkCurvature.setText(_translate("AntennaIntervisibility", "Use earth curvature", None))
        self.txtRefraction.setText(_translate("AntennaIntervisibility", "0.13", None))
        self.label_17.setText(_translate("AntennaIntervisibility", "Atmospheric refraction", None))
        self.chkCumulative.setText(_translate("AntennaIntervisibility", " cumulative (for raster output)", None))
        self.label_7.setText(_translate("AntennaIntervisibility", "Target points (intervisibility)", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.General_tab), _translate("AntennaIntervisibility", "General", None))
        self.textBrowser_2.setHtml(_translate("AntennaIntervisibility", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Ubuntu\'; font-size:10pt; font-weight:600;\">Raster layer:</span><span style=\" font-family:\'Ubuntu\'; font-size:10pt;\"> any supported raster format. For better performance the extent of the raster should be cropped to the analysed area. Too large of a raster will saturate memory.</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Ubuntu\'; font-size:10pt; font-weight:600;\">Observer points:</span><span style=\" font-family:\'Ubuntu\'; font-size:10pt;\"> shapefile containing observer points. The coordinate reference systems of the elevation raster and the observer/target point(s) </span><span style=\" font-family:\'Ubuntu\'; font-size:10pt; text-decoration: underline;\">must match</span><span style=\" font-family:\'Ubuntu\'; font-size:10pt;\">. If field named</span><span style=\" font-family:\'Ubuntu\'; font-size:10pt; font-weight:600;\"> ID </span><span style=\" font-family:\'Ubuntu\'; font-size:10pt;\">exists</span><span style=\" font-family:\'Ubuntu\'; font-size:10pt; font-weight:600;\">, </span><span style=\" font-family:\'Ubuntu\'; font-size:10pt;\">it will be used for filenames and </span><span style=\" font-family:\'Ubuntu\'; font-size:10pt; text-decoration: underline;\">should be unique</span><span style=\" font-family:\'Ubuntu\'; font-size:10pt;\">. Otherwise, the internal Id is used.</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Ubuntu\'; font-size:10pt; font-weight:600;\">Target points:</span><span style=\" font-family:\'Ubuntu\'; font-size:10pt;\"> this option is used for intervisibility analysis only.</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Ubuntu\'; font-size:10pt; font-weight:600;\">Search radius</span><span style=\" font-family:\'Ubuntu\'; font-size:10pt;\">: size of the analyzed area around each observer point.</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Ubuntu\'; font-size:10pt; font-weight:600;\">Adapt to highest point: </span><span style=\" font-family:\'Ubuntu\'; font-size:10pt;\">search the highest point in the vicinity. The search is made in a quadrangular window where the observer point is in the middle.</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Ubuntu\'; font-size:10pt; font-weight:600;\">Viewshed:</span><span style=\" font-family:\'Ubuntu\'; font-size:10pt;\"> standard true/false (binary) viewshed. Multiple viewsheds can be combined in a one raster layer using the cumulative option.</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Ubuntu\'; font-size:10pt; font-weight:600;\">Intervisibility: </span><span style=\" font-family:\'Ubuntu\'; font-size:10pt;\">produces a network of visual relationships between two sets of points (or whithin a single set). </span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Ubuntu\'; font-size:10pt; font-weight:600;\">Invisibility depth:</span><span style=\" font-family:\'Ubuntu\'; font-size:10pt;\"> measures the size an object should attain in order to become visible if placed in an area out of view. </span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Ubuntu\'; font-size:10pt; font-weight:600;\">Horizon:</span><span style=\" font-family:\'Ubuntu\'; font-size:10pt;\"> is the last visible place on the terrain, which corresponds to fringes of visible zones.</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Ubuntu\'; font-size:10pt; font-weight:600;\">Earth curvature:</span><span style=\" font-family:\'Ubuntu\'; font-size:10pt;\"> takes into account the slope of  Earth\'s surface around the observer point.</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Ubuntu\'; font-size:10pt; font-weight:600;\">Atmospheric refraction:</span><span style=\" font-family:\'Ubuntu\'; font-size:10pt;\">  coefficient used to calculate the bending down of the light due to the atmosphere.</span></p></body></html>", None))
        self.label_19.setText(_translate("AntennaIntervisibility", "<html><head/><body><p><a href=\"http://zoran-cuckovic.from.hr/landscape/viewshed-analysis\"><span style=\" text-decoration: underline; color:#0000ff;\">See project web page for more information.</span></a></p></body></html>", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("AntennaIntervisibility", "Reference", None))
        self.textBrowser.setHtml(_translate("AntennaIntervisibility", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:\'Ubuntu\'; font-size:11pt;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Ubuntu\'; font-size:12pt; font-weight:600;\">Version 0.5.1</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:\'Ubuntu\'; font-size:12pt; font-weight:600;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Ubuntu\'; font-size:11pt;\">Developed by Zoran Čučković, Laboratoire Chrono-environnement – UMR 6249, Université de Franche-Comté, Besançon (France).</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:\'Ubuntu\'; font-size:11pt;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Ubuntu\'; font-size:11pt;\">This tool has been brought to you by an individual and is being developed through academic research: please provide an appropriate reference.</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:\'Ubuntu\'; font-size:11pt;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Ubuntu\'; font-size:11pt;\">Contact: cuckovic.zoran@gmail.com</span></p></body></html>", None))
        self.label_20.setText(_translate("AntennaIntervisibility", "<html><head/><body><p><a href=\"http://hub.qgis.org/projects/viewshed/wiki\"><span style=\" text-decoration: underline; color:#0000ff;\">See project web page for more information.</span></a></p></body></html>", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.About_tab), _translate("AntennaIntervisibility", "About", None))

