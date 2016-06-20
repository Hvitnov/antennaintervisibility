# -*- coding: utf-8 -*-
"""
/***************************************************************************
 AntennaIntervisibility
                                 A QGIS plugin
 ------description-------
                              -------------------
        begin                : 2013-05-22
        copyright            : (C) 2013 by Zoran Čučković
        email                : /
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
# Import the PyQt and QGIS libraries

from PyQt4.QtCore import *
from PyQt4.QtGui import *
from qgis.core import *
# Initialize Qt resources from file resources.py
import resources_rc
# Import the code for the dialog
from antennaintervisibilitydialog import AntennaIntervisibilityDialog
from doViewshed import *
import os

    
class AntennaIntervisibility(QObject):
    def __init__(self, iface):
        super(AntennaIntervisibility, self).__init__()
        # Save reference to the QGIS interface
        self.iface = iface
        # initialize plugin directory
        self.plugin_dir = QFileInfo(QgsApplication.qgisUserDbFilePath()).path() + "/python/plugins/antennaintervisibility"
        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'AntennaIntervisibility_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)

            if qVersion() > '4.3.3':
                QCoreApplication.installTranslator(self.translator)

        # Create the dialog (after translation) and keep reference

        # Create the dialog (after translation) and keep reference
        self.dlg = AntennaIntervisibilityDialog()
        self.timeMeasurements = {}

    def initGui(self):
        # Create action that will start plugin configuration
        # icon in the plugin reloader : from resouces.qrc file (compiled)
        self.action = QAction(
            QIcon(":/plugins/AntennaIntervisibility/icon.png"),
            u"Antenna Intervisibility", self.iface.mainWindow())

        # connect the action to the run method
        self.action.triggered.connect(self.run)

        # Add toolbar button and menu item
        self.iface.addToolBarIcon(self.action)
        self.iface.addPluginToMenu(u"&Antenna Intervisibility", self.action)

        # Fire refreshing of combo-boxes containing lists of table columns, after a new layer has been selected
        self.dlg.ui.cmbPoints.currentIndexChanged.connect(self.populate_height_field_comboboxes)
        self.dlg.ui.cmbPointsTarget.currentIndexChanged.connect(self.populate_height_field_comboboxes)

    def unload(self):
        # Remove the plugin menu item and icon
        self.iface.removePluginMenu(u"&Antenna Intervisibility", self.action)
        self.iface.removeToolBarIcon(self.action)

    def populate_height_field_comboboxes(self):
        layer_combobox = self.sender().objectName()
        combobox_index = self.sender().currentIndex()
        layer = str(self.sender().itemData(combobox_index))

        if layer_combobox == 'cmbPoints':
            fields_combobox = self.dlg.ui.cmbObsField
        else:
            fields_combobox = self.dlg.ui.cmbTargetField
        fields_combobox.clear()
        fields_combobox.addItem('', 0)

        map_layer = QgsMapLayerRegistry.instance().mapLayer(layer)

        if map_layer is None: return
        
        provider = map_layer.dataProvider()
        fields = provider.fields() # a dictionary
        for field in fields:
            fields_combobox.addItem(str(field.name()),str(field.name()))

    def printMsg(self, msg):
        QMessageBox.information(self.iface.mainWindow(), "Debug", msg)

    def debugHere(self):
        import pdb
        # These lines allow you to set a breakpoint in the app
        pyqtRemoveInputHook()
        pdb.set_trace()
        return


    def run(self):
        self.timeMeasurements['1 - startTime'] = time.strftime('%H:%M:%S', time.localtime(time.time()))

        #UBACIVANJE RASTERA I TOCAKA (mora biti ovdje ili se barem pozvati odavde)
        myLayers = []
        iface = self.iface
        #clear combos        
        self.dlg.ui.cmbRaster.clear();self.dlg.ui.cmbPoints.clear();self.dlg.ui.cmbPointsTarget.clear()
        #add an empty value to optional combo
        self.dlg.ui.cmbPointsTarget.addItem('',0)
        #add layers to combos
        for i in range(len(iface.mapCanvas().layers())):
            myLayer = iface.mapCanvas().layer(i)
            if myLayer.type() == myLayer.RasterLayer:

                #provjera da li je DEM 1 band .... !!!
                self.dlg.ui.cmbRaster.addItem(myLayer.name(),myLayer.id())

            elif myLayer.geometryType() == QGis.Point: 
                self.dlg.ui.cmbPoints.addItem(myLayer.name(),myLayer.id())
                self.dlg.ui.cmbPointsTarget.addItem(myLayer.name(),myLayer.id())

        #allAttrs = layer.pendingAllAttributesList()
       
                
        # show the dialog
        self.dlg.show()
        # Run the dialog event loop
        result = self.dlg.exec_()

        # See if OK was pressed
        
        if result == 1:
            self.timeMeasurements['2 - oktime'] = time.strftime('%H:%M:%S', time.localtime(time.time()))

            outPath = AntennaIntervisibilityDialog.get_output_file_path(self.dlg)
            raster_layer = AntennaIntervisibilityDialog.get_selected_raster_layer(self.dlg)
            observers_layer = AntennaIntervisibilityDialog.get_selected_observer_layer(self.dlg)
            targets_layer = AntennaIntervisibilityDialog.get_selected_target_layer(self.dlg)

            observer_height = AntennaIntervisibilityDialog.get_observer_height_preset(self.dlg)
            observer_height_field = self.dlg.ui.cmbObsField.itemData(
                self.dlg.ui.cmbObsField.currentIndex())

            target_height = AntennaIntervisibilityDialog.get_target_height_preset(self.dlg)
            target_height_field =self.dlg.ui.cmbTargetField.itemData(
                self.dlg.ui.cmbTargetField.currentIndex())       
            
            Radius = AntennaIntervisibilityDialog.get_observer_view_radius(self.dlg)
            
            output_options = AntennaIntervisibilityDialog.get_algorithm_type(self.dlg)

            
            if not output_options [0]:
                QMessageBox.information(self.iface.mainWindow(), "Error!", str("Select an output option")) 
                return
            self.timeMeasurements['3 - preViewshedTime'] = time.strftime('%H:%M:%S', time.localtime(time.time()))

            out_raster = Viewshed(observers_layer, raster_layer, observer_height, target_height, Radius,outPath,
                                  output_options, self.timeMeasurements,
                                  targets_layer, observer_height_field, target_height_field)

            self.timeMeasurements['4 - preRasterTime'] = time.strftime('%H:%M:%S', time.localtime(time.time()))

            for r in out_raster:
                #QMessageBox.information(self.iface.mainWindow(), "debug", str(r))
                lyName = os.path.splitext(os.path.basename(r))
                layer = QgsRasterLayer(r, lyName[0])
                #if error -> it's shapefile, skip rendering...
                if not layer.isValid():
                    layer= QgsVectorLayer(r,lyName[0],"ogr")
                    
                else:
##                    #rlayer.setColorShadingAlgorithm(QgsRasterLayer.UndefinedShader)
##
##                    #from linfinity.com
##                    extentMin, extentMax = layer.computeMinimumMaximumFromLastExtent( band )
##
##                    # For greyscale layers there is only ever one band
##                    band = layer.bandNumber( layer.grayBandName() ) # base 1 counting in gdal
##                    # We don't want to create a lookup table
##                    generateLookupTableFlag = False
##                    # set the layer min value for this band
##                    layer.setMinimumValue( band, extentMin, generateLookupTableFlag )
##                    # set the layer max value for this band
##                    layer.setMaximumValue( band, extentMax, generateLookupTableFlag )
##
##                    # let the layer know that the min max are user defined
##                    layer.setUserDefinedGrayMinimumMaximum( True )
##
##                    # make sure the layer is redrawn
##                    layer.triggerRepaint()

                    #NOT WORKING 
                    
##                    x = QgsRasterTransparency.TransparentSingleValuePixel()
##                    x.pixelValue = 0
##                    x.transparencyPercent = 100
##                    layer.setTransparentSingleValuePixelList( [ x ] )
                    
                    layer.setContrastEnhancement(QgsContrastEnhancement.StretchToMinimumMaximum)

                    #rlayer.setDrawingStyle(QgsRasterLayer.SingleBandPseudoColor)
                    #rlayer.setColorShadingAlgorithm(QgsRasterLayer.PseudoColorShader)
                    #rlayer.setContrastEnhancementAlgorithm(QgsContrastEnhancement.StretchToMinimumMaximum, False)
                    #rlayer.setTransparency(200)
                    #rlayer.setNoDataValue(0.0)
                self.debugHere()
                QgsMapLayerRegistry.instance().addMapLayer(layer)

            
#    adding csv files ... an attempt
##                    url = QUrl.fromLocalFile(r)
##                    url.addQueryItem('delimiter',',')
##                    url.addQueryItem('xField  <-or-> yField','longitude')
##                    url.addQueryItem('crs','epsg:4723')
##                    url.addQueryItem('wktField','WKT')
## -> Problem
##                    #layer_uri=Qstring.fromAscii(url.toEncoded())
##                    layer_uri= str(url)
##                    layer=QgsVectorLayer(r, lyName[0],"delimitedtext")

#                     QMessageBox.information(None, "File created!", str("Please load file manually (as comma delilmited text)."))





            
