# -*- coding: utf-8 -*-
"""
/***************************************************************************
AntennaIntervisibility
modified Bresenham Algorithm
Written By Jakob Hvitnov
***************************************************************************/
"""

import numpy as np

def bresenham_3d_line_of_sight(observers, targets, raster, obs_height_field,
                               tar_height_field, radius, raster_crs, fresnel=False):
    """bresenham line of sight algorithm"""
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
        except ValueError as e:
            print e

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
        except AttributeError as e:
            print e

        try:
            obs_height = observer.items()[obs_height_field]
            if obs_height is None:
                obs_height = 1.6 # set observer height to person height
            z = obs_height + raster_array[new_obs_y, new_obs_x]
        except(IndexError, TypeError) as e:
            print e

        start = (new_obs_y, new_obs_x, z)

        writePoints.append([(mask_x, mask_y),
                       (mask_x, mask_y + (mask_height * pixelHeight)),
                       (mask_x + (mask_width * pixelWidth) , mask_y),
                       (mask_x + (mask_width * pixelWidth), mask_y + (mask_height * pixelHeight))])

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
