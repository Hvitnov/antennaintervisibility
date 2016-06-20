# -*- coding: utf-8 -*-
"""
/***************************************************************************
 AntennaIntervisibilityDialog
                                 A QGIS plugin
 ------description-------
                             -------------------
        begin                : 2013-05-22
        copyright            : (C) 2013 by Zoran Čučković
        email                : ----
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

from PyQt4 import QtCore, QtGui
from ui_antennaintervisibility import Ui_AntennaIntervisibility
# novo
from qgis.core import * #treba radi QGis.point (srediti ???)
#from qgis.gui import *

import qgis #jer referencira qgis.utils.iface (zapravo treba samo utils ?)
import subprocess, os, sys

class AntennaIntervisibilityDialog(QtGui.QDialog):
   
    def __init__(self):
        
        QtGui.QDialog.__init__(self)
        
        # Set up the user interface from Designer.
        self.ui = Ui_AntennaIntervisibility()
        self.ui.setupUi(self)
        #self.loadLayers()  IZBRISANO
        
        #connections
        self.ui.cmdBrowse.clicked.connect(self.fileOutput)
##        self.ui.cmdAbout.clicked.connect(self.OpenPDFfile)

    def get_selected_observer_layer(self):
        #ovo mu daje varant sa svim podaccima
        l=self.ui.cmbPoints.itemData(self.ui.cmbPoints.currentIndex())
        return str(l)
    
    def get_selected_raster_layer(self):
        #ovo mu daje varant sa svim podaccima
        l=self.ui.cmbRaster.itemData(self.ui.cmbRaster.currentIndex())
        return str(l)
    
    def get_output_file_path(self):
        #ovo mu daje varant sa svim podaccima
        l=self.ui.txtOutput.toPlainText() #inace text()..
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
        return opt

    def fileOutput(self): #problem je ekstenzija!!!!
        homedir = os.path.expanduser('~') # ...works on at least windows and linux. 
        
        filename = QtGui.QFileDialog.getSaveFileName(self, 'Save File', homedir, '*')
        try :
            fname = open(filename, 'w')
            self.ui.txtOutput.clear()
            self.ui.txtOutput.insertPlainText(filename)
            #fname.write(txtOutput.toPlainText())
            fname.close()
        except: pass

    def setProgressBar(self, total):
        self.ui.progressBar.setMinimum(1)
        self.ui.progressBar.setMaximum(total)

    def updateProgressBar(self, val):
        self.ui.progressBar.setValue(val)

    
        
