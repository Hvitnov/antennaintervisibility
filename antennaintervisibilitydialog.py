# -*- coding: utf-8 -*-
"""
/***************************************************************************
AntennaIntervisibility
QGIS plugin
Boilerplate code
***************************************************************************/
"""

from PyQt4 import QtGui
from ui_antennaintervisibility import Ui_AntennaIntervisibility
import os

class AntennaIntervisibilityDialog(QtGui.QDialog):
   
    def __init__(self):
        
        QtGui.QDialog.__init__(self)
        
        # Set up the user interface from Designer.
        self.ui = Ui_AntennaIntervisibility()
        self.ui.setupUi(self)

        self.ui.cmdBrowse.clicked.connect(self.fileOutput)

    def get_selected_observer_layer(self):
        l=self.ui.cmbPoints.itemData(self.ui.cmbPoints.currentIndex())
        return str(l)
    
    def get_selected_raster_layer(self):
        l=self.ui.cmbRaster.itemData(self.ui.cmbRaster.currentIndex())
        return str(l)
    
    def get_output_file_path(self):
        l=self.ui.txtOutput.toPlainText()
        return str(l)
    
    def get_selected_target_layer(self):
        k=self.ui.cmbPointsTarget.currentIndex()
        if k>0:
            l=self.ui.cmbPointsTarget.itemData(k)
            return str(l)
        else: return 0

    def get_observer_height_preset(self):
        try: l = float(self.ui.txtObserver.text())
        except: l=0
        return l

    def get_target_height_preset(self):
        try: l = float(self.ui.txtTarget.text())
        except: l = 0
        return l
    
    def get_observer_view_radius(self):
        try: l = float(self.ui.txtRadius.text())
        except: l = 0
        return l
    
    def get_algorithm_type(self):
        opt = [0,0,0]
        if self.ui.chkIntervisibility.isChecked():  opt [0]= "Intervisibility"
        if self.ui.chkIntervisibility_2.isChecked():  opt[0] = "Intervisibility 2"
        if self.ui.chkIntervisibility_3.isChecked():  opt[0] = "Intervisibility 3"
        return opt

    def fileOutput(self):
        homedir = os.path.expanduser('~')
        filename = QtGui.QFileDialog.getSaveFileName(self, 'Save File', homedir, '*')
        try :
            fname = open(filename, 'w')
            self.ui.txtOutput.clear()
            self.ui.txtOutput.insertPlainText(filename)
            fname.close()
        except: pass

    def setProgressBar(self, total):
        self.ui.progressBar.setMinimum(1)
        self.ui.progressBar.setMaximum(total)

    def updateProgressBar(self, val):
        self.ui.progressBar.setValue(val)

    
        
