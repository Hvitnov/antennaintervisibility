# -*- coding: utf-8 -*-
"""
/***************************************************************************
AntennaIntervisibility
Boilerplate code
QGIS plugin
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
        iface = self.iface
        #clear combos        
        self.dlg.ui.cmbRaster.clear();self.dlg.ui.cmbPoints.clear();self.dlg.ui.cmbPointsTarget.clear()
        #add an empty value to optional combo
        self.dlg.ui.cmbPointsTarget.addItem('',0)
        #add layers to combos
        for i in range(len(iface.mapCanvas().layers())):
            myLayer = iface.mapCanvas().layer(i)
            if myLayer.type() == myLayer.RasterLayer:
                self.dlg.ui.cmbRaster.addItem(myLayer.name(),myLayer.id())

            elif myLayer.geometryType() == QGis.Point: 
                self.dlg.ui.cmbPoints.addItem(myLayer.name(),myLayer.id())
                self.dlg.ui.cmbPointsTarget.addItem(myLayer.name(),myLayer.id())

        # show the dialog
        self.dlg.show()
        # Run the dialog event loop
        result = self.dlg.exec_()

        # See if OK was pressed
        if result == 1:
            self.timeMeasurements['1 - oktime'] = time.strftime('%H:%M:%S', time.localtime(time.time()))
            self.timeMeasurements['1 - oktime_ms'] = time.time()

            outpath = AntennaIntervisibilityDialog.get_output_file_path(self.dlg)
            raster_layer_selection = AntennaIntervisibilityDialog.get_selected_raster_layer(self.dlg)
            observers_layer_selection = AntennaIntervisibilityDialog.get_selected_observer_layer(self.dlg)
            targets_layer_selection = AntennaIntervisibilityDialog.get_selected_target_layer(self.dlg)

            observer_height = AntennaIntervisibilityDialog.get_observer_height_preset(self.dlg)
            observers_height_field = self.dlg.ui.cmbObsField.itemData(self.dlg.ui.cmbObsField.currentIndex())

            target_height = AntennaIntervisibilityDialog.get_target_height_preset(self.dlg)
            targets_height_field = self.dlg.ui.cmbTargetField.itemData(self.dlg.ui.cmbTargetField.currentIndex())
            
            observer_radius = AntennaIntervisibilityDialog.get_observer_view_radius(self.dlg)
            
            output_options = AntennaIntervisibilityDialog.get_algorithm_type(self.dlg)

            raster_map_layer = QgsMapLayerRegistry.instance().mapLayer(raster_layer_selection)
            raster_crs = raster_map_layer.crs()
            raster_path = raster_map_layer.dataProvider().dataSourceUri()
            raster = gdal.Open(raster_path)

            observers_map_layer = QgsMapLayerRegistry.instance().mapLayer(observers_layer_selection)
            observers_path = observers_map_layer.dataProvider().dataSourceUri()
            targets_map_layer = QgsMapLayerRegistry.instance().mapLayer(targets_layer_selection)
            targets_path = targets_map_layer.dataProvider().dataSourceUri()
            if '|' in observers_path:
                path_end = observers_path.find('|')
                observers_path = observers_path[:path_end]
            if '|' in targets_path:
                path_end = targets_path.find('|')
                targets_path = targets_path[:path_end]

            observers_ogr = ogr.Open(observers_path)
            observers = observers_ogr.GetLayer()
            targets_ogr = ogr.Open(targets_path)
            targets = targets_ogr.GetLayer()



            if not output_options [0]:
                QMessageBox.information(self.iface.mainWindow(), "Error!", str("Select an output option")) 
                return

            if self.dlg.ui.chkIntervisibility.isChecked():
                out_raster = Viewshed(observers_layer_selection, raster_layer_selection, observer_height,
                                      target_height, observer_radius,outpath, output_options, targets_layer_selection,
                                      observers_height_field, targets_height_field)

                for r in out_raster:
                    lyName = os.path.splitext(os.path.basename(r))
                    layer = QgsRasterLayer(r, lyName[0])

                    if not layer.isValid():
                        layer = QgsVectorLayer(r, lyName[0], "ogr")
                    else:
                        layer.setContrastEnhancement(QgsContrastEnhancement.StretchToMinimumMaximum)

                    QgsMapLayerRegistry.instance().addMapLayer(layer)

            elif self.dlg.ui.chkIntervisibility_2.isChecked():
                fresnel = self.dlg.ui.checkBox_Fresnel.isChecked()
                bresresult = bresenham_3d_line_of_sight(observers,
                                                        targets,
                                                        raster,
                                                        observers_height_field,
                                                        targets_height_field,
                                                        observer_radius,
                                                        raster_crs,
                                                        fresnel=fresnel)

                # Write line to shapefile
                shpfile = self.write_lines_layer(outpath, raster_map_layer.crs(), bresresult)

                # add shapefile as QGIS layer
                path = os.path.abspath(shpfile)
                basename = os.path.splitext(os.path.basename(shpfile))[0]
                layer = QgsVectorLayer(path, basename, "ogr")
                QgsMapLayerRegistry.instance().addMapLayer(layer)


            self.timeMeasurements['2 - endtime'] = time.strftime('%H:%M:%S', time.localtime(time.time()))
            self.timeMeasurements['2 - endtime_ms'] = time.time()
            elapsedTime = self.timeMeasurements['2 - endtime_ms'] - self.timeMeasurements['1 - oktime_ms']
            self.timeMeasurements['result'] = time.strftime('%H:%M:%S', time.localtime(elapsedTime))
            self.timeMeasurements['result_ms'] = elapsedTime
            self.printMsg(str(self.timeMeasurements))
            return True


    def write_lines_layer(self, file_name, coordinate_ref_system, data_list):
        fields = QgsFields()

        fields.append(QgsField("observ_id", QVariant.String, 'string', 10))
        fields.append(QgsField("target_id", QVariant.String, 'string', 10))
        fields.append(QgsField("visible", QVariant.String, 'string', 10))

        writer = QgsVectorFileWriter(file_name + ".shp", "CP1250", fields,
                                     QGis.WKBLineString, coordinate_ref_system)  # , "ESRI Shapefile"
        # CP... = encoding
        if writer.hasError() != QgsVectorFileWriter.NoError:
            QMessageBox.information(None, "ERROR!", "Cannot write point file - did you select a path?")
            return 0

        for data in data_list:
            # create a new feature
            feat = QgsFeature()
            feat.setFields(fields)
            # Write point data
            start_x  = data['observer_coordinates'][0]
            start_y = data['observer_coordinates'][1]
            end_x, end_y = data['target_coordinates']
            start_point = QgsPoint(start_x, start_y)
            end_point = QgsPoint(end_x, end_y)
            feat.setGeometry(QgsGeometry().fromPolyline([start_point, end_point]))
            feat['observ_id'] = str(data['observer_id'])
            feat['target_id'] = str(data['target_id'])
            feat['visible'] = str(data['visible'])
            writer.addFeature(feat)
            del feat
        del writer
        return file_name + ".shp"
