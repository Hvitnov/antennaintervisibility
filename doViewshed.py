# -*- coding: utf-8 -*-

"""
/***************************************************************************
AntennaIntervisibility
A QGIS plugin
begin : 2013-05-22
copyright : (C) 2013 by Zoran Čučković
email : /
***************************************************************************/

/***************************************************************************
* *
* This program is free software; you can redistribute it and/or modify *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation; either version 2 of the License, or *
* (at your option) any later version. *
* *
***************************************************************************/
"""

#open console !!!! to monitor !!!
#and test time !!!!






from __future__ import division

import qgis
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from qgis.core import *
from qgis.utils import *
from osgeo import osr, gdal, ogr
import time
import numpy
import numpy as np


def debugHere():
    import pdb
    # These lines allow you to set a breakpoint in the app
    pyqtRemoveInputHook()
    pdb.set_trace()
    return

def dist(x1,y1,x2,y2, estimation=False):
    if not estimation: r = numpy.sqrt(numpy.power((x1-x2),2) + numpy.power((y1-y2),2))
    else: # error = cca 1% - NOT USED!
        rt= 1.4142135623730950488016887242097
        r = (numpy.abs(x1-x2) * rt + numpy.abs(y1-y2) * rt) / 2
    return r

def sum_by_group(values, groups):
    #http://stackoverflow.com/questions/4373631/sum-array-by-number-in-numpy
    order = numpy.argsort(groups)
    groups = groups[order]
    values = values[order]
    values.cumsum(out=values)
    index = numpy.ones(len(groups), 'bool')
    index[:-1] = groups[1:] != groups[:-1]
    values = values[index]
    groups = groups[index]
    values[1:] = values[1:] - values[:-1]
    return values, groups


def get_fresnel_radius(distance_to_observer, distance_to_target):
    speedOfLight = 300000000.0
    frequency_wifi = 2400000000.0
    wavelength = speedOfLight / frequency_wifi

    # distance from obs
    d1 = distance_to_observer
    # distance from target
    d2 = distance_to_target
    # Fresnel xone number
    zone = 1

    # Fresnel radius formula
    fresnel_radius = np.sqrt(
        (zone*wavelength*d1*d2)/
        (d1 + d2))

    return fresnel_radius

def euclidian_distance(x,y):
    if type(x) == tuple:
        x = np.array(x)
    if type(y) == tuple:
        y = np.array(y)
    return np.sqrt(np.sum((x-y)**2))


def bresenham_3d_line_of_sight(observers, targets, raster, obs_height_field,
                               tar_height_field, radius, raster_crs, fresnel=False):
    """Naive bresenham line of sight algorithm"""
    writePoints = []
    lines_for_shp = []
    start, end = 0, 0

    raster_transform = raster.GetGeoTransform()
    pixelWidth = raster_transform[1]
    pix = pixelWidth
    pixelHeight = raster_transform[5]
    xOrigin = raster_transform[0]
    yOrigin = raster_transform[3]

    raster_band = raster.GetRasterBand(1)
    info = []
    for obs in range(observers.GetFeatureCount()):
        observer = observers.GetFeature(obs)
        # get Observer point geometry
        obs_geom = observer.geometry()
        try:
            obs_x = obs_geom.GetPoints()[0][0]
            obs_y = obs_geom.GetPoints()[0][1]
        except ValueError:
            debugHere()
        # offset x,y values to equivalent raster index values
        obs_x_off = int((obs_x - xOrigin) / pixelWidth)
        obs_y_off = int((obs_y - yOrigin) / pixelHeight)
        mask_x = obs_x - radius
        mask_y = obs_y - radius
        mask_x_pix = int((mask_x - xOrigin) / pixelWidth)
        mask_y_pix = int((mask_y - yOrigin) / pixelHeight)
        radius_pix = int(radius / pixelWidth)
        mask_width = radius_pix * 2
        mask_height = radius_pix * 2
        if mask_x_pix < 0: # mask has overflow beyond raster edge
            mask_width += mask_x_pix # clip mask width by the overflow
            mask_x_pix = 0 # set mask origin x to edge of raster
            mask_x = xOrigin
        if mask_y_pix < 0:
            mask_height += mask_y_pix
            mask_y_pix = 0
            mask_y = yOrigin
        # truncate positive overflow
        if mask_width + mask_x_pix > raster_band.XSize:
            overflow = raster_band.XSize - (mask_width + mask_x_pix)
            mask_width += overflow
        if mask_height + mask_y_pix > raster_band.YSize:
            overflow = raster_band.YSize - (mask_height + mask_y_pix)
            mask_height += overflow
        mask_x_pix = int(mask_x_pix)
        mask_y_pix = int(mask_y_pix)
        mask_width = int(mask_width)
        mask_height = int(mask_height)
        new_obs_x = obs_x_off - mask_x_pix
        new_obs_y = mask_y_pix - obs_y_off

        # x_geog, y_geog = raster_x_min + x * pix + pix / 2, raster_y_max - y * pix - pix / 2
        # areaOfInterest = QgsRectangle(x_geog - radius, y_geog - radius, x_geog + radius, y_geog + radius)
        # set observer height

        # Raster used is smaller than radius, so no clipping nescesarry
        try:
            if raster_band.YSize < radius * 2 or raster_band.YSize < radius * 2:
                mask_x = xOrigin
                mask_y = yOrigin
                new_obs_x = obs_x_off
                new_obs_y = obs_y_off
                raster_array = raster_band.ReadAsArray().astype(np.float)
            else:
                raster_array = raster_band.ReadAsArray(mask_x_pix, mask_y_pix, mask_width, mask_height).astype(np.float)
        except:
            debugHere()
        try:
            obs_height = observer.items()[obs_height_field]
            if obs_height is None:
                obs_height = 1.6 # set observer height to person height
            z = obs_height + raster_array[new_obs_y, new_obs_x]
        except(IndexError, TypeError) as e:
            print e
            debugHere()
        start = (new_obs_y, new_obs_x, z)

        writePoints.append([(mask_x, mask_y),
                       (mask_x, mask_y + (mask_height * pixelHeight)),
                       (mask_x + (mask_width * pixelWidth) , mask_y),
                       (mask_x + (mask_width * pixelWidth), mask_y + (mask_height * pixelHeight))])

        # raster_crs

        for tar in range(targets.GetFeatureCount()):
            target_in_radius = True
            target = targets.GetFeature(tar)
            # get Target point geometry
            tar_geom = target.geometry()
            x, y = tar_geom.GetPoints()[0]

            target_outside_radius = euclidian_distance((obs_x, obs_y), (x, y)) > radius
            if target_outside_radius:
                continue
            # offset x,y values to equivalent raster index values
            x = int((x - mask_x) / pixelWidth)
            y = int((y - mask_y) / pixelHeight)




            # check if target point is out of search area
            # if target_outside_radius:
            #    continue

            # get target height
            z = target.items()[tar_height_field]

            try:
                landscape_height = raster_array[y, x]
            except IndexError:
                target_in_radius = False
                continue

            # get target height
            z = target.items()[tar_height_field] + landscape_height
            end = (y, x, z)

            # Unpack start/end tuples
            x, y, z = start
            x2, y2, z2 = end
            z_value = z

            # Calculate differentials
            diff_x = x2 - x
            diff_y = y2 - y
            diff_z = z2 - z

            # Assign incremental slope values for x, y, z
            incr_x = -1 if (diff_x < 0) else 1
            incr_y = -1 if (diff_y < 0) else 1
            incr_z = -1 if (diff_z < 0) else 1

            abs_diff_x = abs(diff_x)
            abs_diff_y = abs(diff_y)
            abs_diff_z = abs(diff_z)

            diff_x2 = abs_diff_x * 2
            diff_y2 = abs_diff_y * 2
            diff_z2 = abs_diff_z * 2

            # Find the steepest axis and find line segments accordingly
            if (abs_diff_x >= abs_diff_y) and (abs_diff_x >= abs_diff_z):
                steepest = 'x'
                z_line_length = np.sqrt(pow(diff_x, 2) + pow(diff_z, 2))
                z_segment_length = z_line_length / diff_x
            elif (abs_diff_y > abs_diff_x) and (abs_diff_y >= abs_diff_z):
                steepest = 'y'
                z_line_length = np.sqrt(pow(diff_y, 2) + pow(diff_z, 2))
                z_segment_length = z_line_length / diff_y
            elif (abs_diff_z > abs_diff_x) and (abs_diff_z > abs_diff_y):
                steepest = 'z'
                z_line_length = np.sqrt(pow(diff_x, 2) + pow(diff_z, 2))
                z_segment_length = z_line_length / diff_z
            else:
                return "Error when finding steepest line"

            incr_z_value = np.sqrt(abs(pow(z_segment_length, 2) - pow(1, 2)))
            incr_z_value = -incr_z_value if (diff_z < 0) else incr_z_value

            xm, ym, zm = (x2 + x) / 2, (y2 + y) / 2, (z2 + z) / 2
            zm = z + xm * incr_z_value

            mid_fresnel = get_fresnel_radius(z_line_length / 2, z_line_length / 2)

            if fresnel:
                try:
                    visibility = zm - mid_fresnel > raster_array[xm, ym]
                except:
                    debugHere()
                if not visibility:
                    lines_for_shp.append(build_return_package(observer, target, visibility))
                    continue

            if 'x' in steepest:
                err_1 = diff_y2 - abs_diff_x
                err_2 = diff_z2 - abs_diff_x

                for i in np.arange(abs_diff_x - 1):
                    if (err_1 > 0):
                        y += incr_y
                        err_1 -= diff_x2

                    if (err_2 > 0):
                        z += incr_z
                        err_2 -= diff_x2

                    err_1 += diff_y2
                    err_2 += diff_z2
                    x += incr_x
                    z_value += incr_z_value
                    visibility = z_value > raster_array[x, y]
                    if not visibility:
                        break

            if 'y' in steepest:
                err_1 = diff_x2 - abs_diff_y
                err_2 = diff_z2 - abs_diff_y

                for i in np.arange(abs_diff_y - 1):

                    if (err_1 > 0):
                        x += incr_x
                        err_1 -= diff_y2

                    if (err_2 > 0):
                        z += incr_z
                        err_2 -= diff_y2

                    err_1 += diff_x2
                    err_2 += diff_z2
                    y += incr_y
                    z_value += incr_z_value
                    visibility = z_value > raster_array[x, y]
                    if not visibility:
                        break

            if 'z' in steepest:
                err_1 = diff_y2 - abs_diff_z
                err_2 = diff_x2 - abs_diff_z

                for i in np.arange(abs_diff_z - 1):

                    if (err_1 > 0):
                        y += incr_y
                        err_1 -= diff_z2

                    if (err_2 > 0):
                        x += incr_x
                        err_2 -= diff_z2

                    err_1 += diff_y2
                    err_2 += diff_x2
                    z += incr_z
                    z_value += incr_z_value
                    visibility = z_value > raster_array[x, y]
                    if not visibility:
                        break

            lines_for_shp.append(build_return_package(observer, target, visibility))
    return lines_for_shp

def build_return_package(observer, target, visibility):
    package = {}
    # if features have id fields use them - otherwise use FID
    obs_id_index = observer.GetFieldIndex('ID')
    tar_id_index = target.GetFieldIndex('ID')
    if obs_id_index >= 0:
        package['observer_id'] = observer.GetField(obs_id_index)
    else:
        package['observer_id'] = observer.GetFID()
    if tar_id_index >= 0:
        package['target_id'] = target.GetField(tar_id_index)
    else:
        package['target_id'] = target.GetFID()

    # add coordinates
    package['observer_coordinates'] = observer.geometry().GetPoints()[0]
    package['target_coordinates'] = target.geometry().GetPoints()[0]

    # add visibility
    package['visible'] = visibility
    return package

def writeToPolygonFile(inputPoints, filename, crs):
    pointsList = []
    for pLine in inputPoints:
        polygonPoints = []
        for pt in pLine:
            x, y = pt
            polygonPoints.append(QgsPoint(x, y))
        pointsList.append(polygonPoints)

    fields = QgsFields()
    outfile = "/home/hvitnov/" + filename
    writer = QgsVectorFileWriter(outfile + ".shp", "CP1250", fields, QGis.WKBPoint, crs)

    for points in pointsList:
        feat = QgsFeature()
        feat.setGeometry(QgsGeometry.fromPolygon([points]))
        feat.setFields(fields)
        writer.addFeature(feat)
        del feat

    del writer


def Viewshed (Obs_points_layer, Raster_layer, z_obs, z_target, radius, output,
              output_options,
              Target_layer = 0, z_obs_field = 0, z_target_field = 0): #for visib index !!




# -------------- MAIN ----------------------------

    #####################################
    # Prepare  progress Bar
    
    iface.messageBar().clearWidgets()
    progressMessageBar = iface.messageBar()
    progress_bar = QProgressBar()

    #get a vector layer
    vlayer = QgsMapLayerRegistry.instance().mapLayer(Obs_points_layer)
    vlayer.selectAll()
    #Count all selected feature
    feature_count = vlayer.selectedFeatureCount() 

    #could be set to 100, making it easy to work with percentage of completion
    progress_bar.setMaximum(feature_count) 
    #pass the progress bar to the message Bar
    progressMessageBar.pushWidget(progress_bar)

    #set a counter to reference the progress, 1 will be given after pre-calculations
    progress = 0

    ###########################"
    #read data, check etc
    out_files=[];rpt=[];connection_list=[]

    RasterPath= str(QgsMapLayerRegistry.instance().mapLayer(Raster_layer).dataProvider().dataSourceUri())

    gdal_raster = gdal.Open(RasterPath)
    gt = gdal_raster.GetGeoTransform()
    projection = gdal_raster.GetProjection()
    global pix; pix=gt[1] # there's a bug in Python: globals cannot be avoided by a nested function (??)
    global raster_x_min; raster_x_min = gt[0]
    global raster_y_max; raster_y_max = gt[3] # it's top left y, so maximum!

    raster_y_size = gdal_raster.RasterYSize
    raster_x_size = gdal_raster.RasterXSize

    raster_y_min = raster_y_max - raster_y_size * pix
    raster_x_max = raster_x_min + raster_x_size * pix


    #adfGeoTransform[0] /* top left x */
    #adfGeoTransform[1] /* w-e pixel resolution */
    #adfGeoTransform[2] /* rotation, 0 if image is "north up" */
    #adfGeoTransform[3] /* top left y */
    #adfGeoTransform[4] /* rotation, 0 if image is "north up" */
    #adfGeoTransform[5] /* n-s pixel resolution */

    gtiff = gdal.GetDriverByName('GTiff')#for creting new rasters
    #projection = gdal_raster.GetProjection() #not good ? better to use crs from the project (may override the original)
 #   rb=gdal_raster.GetRasterBand(1)#treba li to?


    Obs_layer=QgsMapLayerRegistry.instance().mapLayer(Obs_points_layer)
    if Obs_layer.isValid():
        # returns 0-? for indexes or -1 if doesn't exist
        obs_has_ID = bool( Obs_layer.fieldNameIndex ("ID") + 1)
    else: return # abandon function       

    #initialise target points and create spatial index for speed
    if Target_layer:
        
        Target_layer = QgsMapLayerRegistry.instance().mapLayer(Target_layer)
        if Target_layer.isValid():
            
            targ_has_ID = bool(Target_layer.fieldNameIndex ("ID") + 1)
          
            targ_index = QgsSpatialIndex()
            for f in Target_layer.getFeatures():
                targ_index.insertFeature(f)
                
        else: return # abandon function
    ################################################
    # precalculate distance matrix, errors etc 
    radius_pix = int(radius/pix)
       
    full_window_size = radius_pix *2 + 1
    
    mx_vis = numpy.zeros((full_window_size, full_window_size))

    temp_x= ((numpy.arange(full_window_size) - radius_pix) * pix) **2
    temp_y= ((numpy.arange(full_window_size) - radius_pix) * pix) **2

    mx_dist = numpy.sqrt(temp_x[:, None] + temp_y[None, :])
    mask_circ = mx_dist [:] > radius  #it's real (metric) radius


    # ----------------- POINT LOOP -------------------------
    for feat in Obs_layer.getFeatures():
        targets=[]

        geom = feat.geometry()
        t = geom.asPoint()
        
        id1 = feat["ID"] if obs_has_ID else feat.id()
        x_geog, y_geog = t[0], t[1]
        
        #check if the point is out of raster extents
        if raster_x_min <= x_geog <= raster_x_max and raster_y_min <= y_geog <= raster_y_max:
            x = int((x_geog - raster_x_min) / pix) # not float !
            y = int((raster_y_max - y_geog) / pix) #reversed !
        else: continue
            
        # find correct coordinates of new points - needed for intervisibilty only...
        if output_options[0]== "Intervisibility":
            x_geog, y_geog= raster_x_min  + x *pix + pix/2 , raster_y_max - y *pix - pix/2
        else: z = 0 #reset !!

        # ----------  extraction of a chunk of data ---------------

        if x <= radius_pix:  #cropping from the front
            x_offset =0
            x_offset_dist_mx = radius_pix - x
        else:               #cropping from the back
            x_offset = x - radius_pix
            x_offset_dist_mx= 0
                            
        x_offset2 = min(x + radius_pix +1, raster_x_size) #could be enormus radius, so check both ends always
            
        if y <= radius_pix:
            y_offset =0
            y_offset_dist_mx= radius_pix - y
        else:
            y_offset = y - radius_pix
            y_offset_dist_mx= 0

        y_offset2 = min(y + radius_pix + 1, raster_y_size )

        window_size_y = y_offset2 - y_offset
        window_size_x = x_offset2 - x_offset

        data = numpy.zeros((full_window_size, full_window_size)) #window is always the same size
        
        data[ y_offset_dist_mx : y_offset_dist_mx +  window_size_y,
              x_offset_dist_mx : x_offset_dist_mx + window_size_x] = gdal_raster.ReadAsArray(
                  x_offset, y_offset, window_size_x, window_size_y).astype(float)# global variable

        if z_obs_field:
            try:    z = data [radius_pix,radius_pix] + float(feat[z_obs_field])
            except: z = data [radius_pix,radius_pix] +  z_obs 
        else:	    z = data [radius_pix,radius_pix] + z_obs  
        
        data -= z # level all according to observer

        # ------  create an array of additional angles (parallel to existing surface) ----------
        if z_target > 0 :
            mx_target = (data + z_target) / mx_dist 

        else: mx_target=None

        data /= mx_dist #all one line = (data -z - mxcurv) /mx_dist
                
        if Target_layer and output_options[0]== "Intervisibility":
          
            areaOfInterest= QgsRectangle (x_geog -radius , y_geog -radius, x_geog +radius, y_geog +radius)
            # SPATIAL INDEX FOR SPEED
            feature_ids = targ_index.intersects(areaOfInterest)

            visib_list=[] #initialize output list
                   
            for fid in feature_ids: # NEXT TO GET THE FEATURE
                feat2 = Target_layer.getFeatures(QgsFeatureRequest().setFilterFid(fid)).next()
                id2 = feat2["ID"] if targ_has_ID else feat2.id()
                geom2 = feat2.geometry()
                x2_geog, y2_geog = geom2.asPoint()
                
                # they may fall out of the entire raster
                if raster_x_min <= x2_geog <= raster_x_max and raster_y_min <= y2_geog <= raster_y_max:
                    #re-align to relative pixel coords of the data matrix
                    x2 = int((x2_geog - raster_x_min) / pix)  #round vraca float!!
                    y2 = int((raster_y_max - y2_geog) / pix)  #pazi: obratno!!
                else: continue

                #skipping 1
                if x2==x and y2==y : continue

                x2_local = radius_pix + (x2 - x)
                y2_local = radius_pix + (y2 - y)#x and y are also global

                try : #skipping 2
                    if mx_dist [y2_local,x2_local] > radius: continue
                except : continue #out of local raster

                tg_offset = z_target
                if z_target_field: #this is a clumsy addition so that each point might have it's own height
                    try: tg_offset = float(feat2[z_target_field])
                    except: pass
                     
                visib, hgt, d, err = intervisibility_zoran(radius_pix, radius_pix, x2_local, y2_local, tg_offset, id2, mx_dist, data)

                 
                connection_list.append([id1, x_geog ,y_geog, id2, x2_geog, y2_geog,
                                        visib, hgt, d])

                
##            elif output_options[0]== "Intervisibility":
##                #to make it fast : for each pair choose matrix in-between
##                #e.g. linspace + stretch or y = m * x[:, np.newaxis] + b (y =m*x +b...)
##                for x2,y2,z2,id2,x2_geo, y2_geo in targets:
##                    
##                    d=mx_dist[y2,x2]
##                    z_angle = z2/d if z2 else 0
##                    vis  = mx_vis[y2,x2]  + z_angle #z2 is only a difference - main target angle is set up as usual
##                    z_object= z_target+ z2
##                    
##                    hgt= vis * d
##                    if hgt > z_object: hgt = z_object
##                          
##                    visib_list.append([id2, x2_geo, y2_geo, vis>=0,  hgt, d])#, err=0

                    
        ######################################
        #Update the progress bar: point loop


        progress += 1
        progress_bar.setValue(progress) #(progress / feature_count) * 100 = percentage - losing time :)	

        start_etape=time.clock()
        
    #####################################
    if output_options[0]=="Intervisibility":
        success = write_intervisibility_line (output, connection_list, Obs_layer.crs())
        if success : out_files.append(success)

        else : QMessageBox.information(None, "Error writing file !", str(output + '_intervisibility cannot be saved'))

    matrix_final = None; data = None; connections_list=None; v=None; vis=None
    
    iface.messageBar().clearWidgets()

    return out_files


def intervisibility_zoran(x, y, x2, y2, target_offset, id_target, mx_dist, data,
                          interpolate=True):  # x0, y0 are not needed!
    d = mx_dist[y2, x2]  # do before swapping

    dx = abs(x2 - x)
    dy = abs(y2 - y)
    steep = (dy > dx)
    # direction of movement : plus or minus in the coord. system
    sx = 1 if (x2 - x) > 0 else -1
    sy = 1 if (y2 - y) > 0 else -1

    if steep:  # if the line is steep: change x and y
        # x,y = y,x they are the same !!
        x2, y2 = y2, x2
        dx, dy = dy, dx
        sx, sy = sy, sx

    D = 0

    # for interpolation
    # slope = dy / dx *sx *sy #!! the operators for negative quadrants (we do need dx, dy as absolute to verify steepness, to search distances...)

    # begins with the minimum possible angle for the observer pix,
    # thus the first pixel next to observer is always visible!
    visib = True
    angle_block = None
    angle_block2 = angle_block
    angle_hor_max = angle_block
    angle_hor_min = angle_block

    for i in numpy.arange(0, dx):
        ##                y_dec = slope * (x-x2) + y2
        ##                y = int(round(y_dec)) #it considers rounded as float...
        ##                err = y - y_dec
        ##
        ##                interpolated = y-1 if err > 0 else y+1 #inverse : if the good pixel is above the line then search down
        ##                                                    # we don't need to y+sy because the slope has already been multiplied by sy
        # ---- Bresenham's algorithm (understandable variant)
        # http://www.cs.helsinki.fi/group/goa/mallinnus/lines/bresenh.html
        x += sx
        if 2 * (D + dy) < dx:
            D += dy  # y remains
        else:
            y += sy
            D += dy - dx

        # unfortunately, numpy allows for negative values...
        # when x or y are too large, the break is on try-except below
        if x < 0 or y < 0: break

        # --------- unpack coordinates ------------
        if not steep:
            x_pix, y_pix = x, y
        else:
            x_pix, y_pix = y, x

        try:
            angle = data[y_pix, x_pix]
        except:
            break

        angle_exact = 0  # initiate/reset

        if D:
            err = numpy.abs(D / dx)  # D can be a very good substitute (2-3% differences)

            sD = -1 if D < 0 else 1  # inverse : if the good pixel is above the line then search down
            interpolated = y + sD * sy

            # if interpolate<0: break

            if not steep:
                x_pix_interp, y_pix_interp = x, interpolated
            else:
                x_pix_interp, y_pix_interp = interpolated, x

            try:
                angle2 = data[y_pix_interp, x_pix_interp]
            except:
                break

        else:
            err = 0
            angle2 = angle

        # it is not correct to test non-interpolated against interpolated horizon angle!!
        # ... nor to test max>max and min>min, because of possible interpolation shift!
        # only : min angle > max old_angle (and vice versa...)

        if angle < angle_hor_min and angle2 < angle_hor_min:
            visib = False
        elif angle > angle_hor_max and angle2 > angle_hor_max:
            visib = True
        else:  # optimisation: interpolate only when really needed

            angle_exact = angle + (angle2 - angle) * err
            if not angle_hor_exact:
                angle_hor_exact = angle_block + (angle_block2 - angle_block) * block_err

            visib = (angle_exact > angle_hor_exact)

        # catch old values
        if visib:
            angle_block, angle_block2, block_err = angle, angle2, err
            angle_hor_exact = angle_exact  # mostly is 0 !

            if angle > angle2:
                angle_hor_max, angle_hor_min = angle, angle2

            else:
                angle_hor_max, angle_hor_min = angle2, angle


            # --------------- processing output ----------------

    if visib:  # there is no ambiguity, visible!
        hgt = target_offset

    else:
        # repeat to calculate height (even without target
        angle_off = target_offset / d if target_offset else 0

        angle += angle_off  # ERROR (err) IS NOT POSSIBLE !
        # here it is probable..
        if not angle_hor_exact:
            angle_hor_exact = angle_block + (angle_block2 - angle_block) * block_err

        hgt = (angle - angle_hor_exact) * d
        visib = (hgt > 0)

        if hgt > target_offset: hgt = target_offset
        # in case when the pixel is preceeded by an invisible one and becomes visible only on exact horizon angle,
        # hgt can get augmented by the relative target pixel height

    return [visib, hgt, d, err]




#works with global gdal_raster
def write_raster (matrix, file_name,columns_no, rows_no , offset_x, offset_y,
                  geotransform_data, GDAL_projection_data,
                  num_format=gdal.GDT_Float32): #full file path

    driver = gdal.GetDriverByName( 'GTiff' )
    
    dst_ds = driver.Create( file_name+'.tiff', columns_no, rows_no, 1, num_format)
    if not dst_ds: return 0
       
    dst_ds.SetProjection(GDAL_projection_data)
    dst_ds.SetGeoTransform(geotransform_data)
    
    dst_ds.GetRasterBand(1).Fill(numpy.nan)#this is for esthetic purpose chiefly 
    dst_ds.GetRasterBand(1).SetNoDataValue(numpy.nan)# nans are set for all outputs if not cumulative
    
    dst_ds.GetRasterBand(1).WriteArray(matrix,offset_x,offset_y)#offset=0
    
    dst_ds=None
    #driver.GDALClose(dst_ds) not working!
    return file_name +'.tiff'

def write_intervisibility_line (file_name, data_list, coordinate_ref_system, use_pix_coords=False):

    #QMessageBox.information(None, "Timing report:", str(data_list))
    
    fields = QgsFields() #there's a BUG in QGIS here (?), normally : fields = ....
    fields.append(QgsField("Source", QVariant.String ))
    fields.append(QgsField("Target", QVariant.String))
## fields.append(QgsField("Source_lbl", QVariant.String, 'string',50))
## fields.append(QgsField("Target_lbl", QVariant.String, 'string',50))
    fields.append(QgsField("Visible", QVariant.String, 'string',5))
    fields.append(QgsField("TargetSize", QVariant.Double, 'double',10,3))
    fields.append(QgsField("Distance", QVariant.Double, 'double',10,2))

    writer = QgsVectorFileWriter( file_name + ".shp", "CP1250", fields,
                                  QGis.WKBLineString, coordinate_ref_system) #, "ESRI Shapefile"
                                            #CP... = encoding
    if writer.hasError() != QgsVectorFileWriter.NoError:
        QMessageBox.information(None, "ERROR!", "Cannot write intervisibilty file (?)")
        return 0
    
    for r in data_list:
        # create a new feature
        feat = QgsFeature()
        if use_pix_coords:
            half_pix= pix/2 #global variable pix
            l_start=QgsPoint(raster_x_min  + r[1]*pix + half_pix, raster_y_max - r[2]*pix - half_pix )
            l_end = QgsPoint(raster_x_min  + r[4]*pix + half_pix, raster_y_max - r[5]*pix - half_pix)
        else:
            l_start=QgsPoint(r[1],r[2]);  l_end = QgsPoint(r[4],r[5])
         
        feat.setGeometry(QgsGeometry.fromPolyline([l_start, l_end]))
        # do not cast ID to string: unicode problem -- angle * distance in pixels -- distance * pixel_size
        #feat.setAttributes([ str(r[0]), str(r[3]), bool(r[6]), float(r[7] * r[8]), ])
        feat.setFields(fields)
        feat['Source'] = r[0]
        feat['Target'] = r[3]
        feat['Visible'] = 'True' if r[6] else 'False'
        feat['TargetSize'] = float(r[7])
        feat['Distance'] = float(r[8])
        
        writer.addFeature(feat)
        del feat

    del writer
    layer = None
    return file_name + ".shp"


