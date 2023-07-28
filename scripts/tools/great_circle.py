# Function to work out a map of angular distances from a pixel

def distancemap(pixel_lon, pixel_lat, nside, good_pixels, lon, lat):
    ''' Return a map of angular distances in degrees from a
    given pixel. All inputs and outputs are in degrees'''

    import numpy as np

    # Initialise distance map as map of np.nans
    distance_map = np.zeros((12*nside**2))*np.nan

    # Parallelized Haversine formula goes here
    def haversine(lon1, lat1, lon2, lat2):
        ''' Calculate the angular distance between two points in the sky '''

        # Convert degrees to radians
        lon1, lat1, lon2, lat2 = map(np.deg2rad, [lon1, lat1, lon2, lat2])

        # Haversine formula
        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
        c = 2 * np.arcsin(np.sqrt(a))
        theta = np.rad2deg(c) # Angular distance in degrees
        return theta

    # Fill in distance map
    distance_map[good_pixels] = np.array(haversine(lon1=pixel_lon, lat1=pixel_lat, lon2=lon[good_pixels], lat2=lat[good_pixels]))

    return distance_map



def great_circle_angles(ptlon1, ptlat1, ptlon2, ptlat2, good_pixels, nside):
    ''' Accepts lon and lat in degrees. "ptlon2" and "ptlat2" can be lists of
        values, but ptlon1 and ptlat1 have to be a single number. '''

    import numpy as np

    # Use longdoubles for precision
    ptlon1 = np.longdouble(ptlon1)
    ptlat1 = np.longdouble(ptlat1)
    ptlon2 = np.longdouble(ptlon2)
    ptlat2 = np.longdouble(ptlat2)

    # Initialise great angle maps
    angles_initial = np.zeros((12*nside**2))*np.nan
    angles_final = np.zeros((12*nside**2))*np.nan

    # Only want good pixels
    ptlon2 = ptlon2[good_pixels]
    ptlat2 = ptlat2[good_pixels]

    # Convert to radians for calculations
    ptlon1_radians = np.deg2rad(ptlon1)
    ptlat1_radians = np.deg2rad(ptlat1)
    ptlon2_radians = np.deg2rad(ptlon2)
    ptlat2_radians = np.deg2rad(ptlat2)

    distance_radians=2*np.arcsin(np.sqrt(np.power((np.sin((ptlat1_radians-ptlat2_radians)/2)),2) + np.cos(ptlat1_radians)*np.cos(ptlat2_radians)*np.power((np.sin((ptlon1_radians-ptlon2_radians)/2)),2)))

    # Right now 1e-9 of a full circle is about 0.001 arcsec
    path_fraction = 1e-9
    fractions_to_compute = [path_fraction,1-path_fraction]

    # Initialise latitudes and longitudes along great circle
    npoints = len(ptlon2_radians)
    mylats = np.zeros((2,npoints))
    mylons = np.zeros((2,npoints))

    for i,f in enumerate(fractions_to_compute):
            # f is expressed as a fraction along the route from point 1 to point 2
            A=np.sin((1-f)*distance_radians)/np.sin(distance_radians)
            B=np.sin(f*distance_radians)/np.sin(distance_radians)
            x = A*np.cos(ptlat1_radians)*np.cos(ptlon1_radians) + B*np.cos(ptlat2_radians)*np.cos(ptlon2_radians)
            y = A*np.cos(ptlat1_radians)*np.sin(ptlon1_radians) +  B*np.cos(ptlat2_radians)*np.sin(ptlon2_radians)
            z = A*np.sin(ptlat1_radians) + B*np.sin(ptlat2_radians)
            newlat=np.arctan2(z,np.sqrt(np.power(x,2)+np.power(y,2)))
            newlon=np.arctan2(y,x)
            newlat_degrees = np.degrees(newlat)
            newlon_degrees = np.degrees(newlon)
            mylats[i,:] = newlat_degrees
            mylons[i,:] = newlon_degrees

    # Calculate initial and final displacements from which the angles can be calculated

    # Initial - need to fold longitude differences so that they are around 0
    differences_initial_lon = np.subtract(mylons[0][:], ptlon1)
    differences_initial_lat = np.subtract(mylats[0][:], ptlat1)
    differences_initial_lon[differences_initial_lon>180] = differences_initial_lon[differences_initial_lon>180]-360.
    differences_initial_lon[differences_initial_lon<-180] = differences_initial_lon[differences_initial_lon<-180]+360.

    # Final - need to fold longitude differences so that they are around 0
    differences_final_lon = np.subtract( ptlon2, mylons[1][:] )
    differences_final_lat = np.subtract( ptlat2, mylats[1][:] )
    differences_final_lon[differences_final_lon>180] = differences_final_lon[differences_final_lon>180]-360.
    differences_final_lon[differences_final_lon<-180] = differences_final_lon[differences_final_lon<-180]+360.

    def angle_from_difference(differences_lon, differences_lat):
        ''' Input and output both in degrees '''
        differences_lon = np.deg2rad(differences_lon)
        differences_lat = np.deg2rad(differences_lat)
        angles = np.arctan2(differences_lat, differences_lon) * 180. / np.pi
        return angles

    # Get angles
    angles_initial[good_pixels] = angle_from_difference(differences_initial_lon, differences_initial_lat)
    angles_final[good_pixels] = angle_from_difference(differences_final_lon, differences_final_lat)

    return angles_initial, angles_final



def greatcircle(ptlon1, ptlat1, ptlon2, ptlat2, distance=0.5):
    ''' Returns a set of points along the great circle, where distance determines
    the distance in degrees between subsequent points. All the input coordinates
    are in degrees, and so is the output. '''

    import numpy as np
    import math

    # Work out the distance between the two original points
    ptlon1_radians = math.radians(ptlon1)
    ptlat1_radians = math.radians(ptlat1)
    ptlon2_radians = math.radians(ptlon2)
    ptlat2_radians = math.radians(ptlat2)
    distance_radians=2*math.asin(math.sqrt(math.pow((math.sin((ptlat1_radians-ptlat2_radians)/2)),2) + math.cos(ptlat1_radians)*math.cos(ptlat2_radians)*math.pow((math.sin((ptlon1_radians-ptlon2_radians)/2)),2)))

    # Work out the number of segments
    numberofsegments = int(round((distance_radians/(2*math.pi)*360)/distance)) # num_of_segments
    onelessthansegments = numberofsegments - 1
    fractionalincrement = (1.0/onelessthansegments)

    # Write the starting coordinates
    mylats = []
    mylons = []
    mylats.append(ptlat1)
    mylons.append(ptlon1)

    f = fractionalincrement
    icounter = 1
    while (icounter <  onelessthansegments):
            icountmin1 = icounter - 1
            mylats.append([])
            mylons.append([])
            # f is expressed as a fraction along the route from point 1 to point 2
            A=math.sin((1-f)*distance_radians)/math.sin(distance_radians)
            B=math.sin(f*distance_radians)/math.sin(distance_radians)
            x = A*math.cos(ptlat1_radians)*math.cos(ptlon1_radians) + B*math.cos(ptlat2_radians)*math.cos(ptlon2_radians)
            y = A*math.cos(ptlat1_radians)*math.sin(ptlon1_radians) +  B*math.cos(ptlat2_radians)*math.sin(ptlon2_radians)
            z = A*math.sin(ptlat1_radians) + B*math.sin(ptlat2_radians)
            newlat=math.atan2(z,math.sqrt(math.pow(x,2)+math.pow(y,2)))
            newlon=math.atan2(y,x)
            newlat_degrees = math.degrees(newlat)
            newlon_degrees = math.degrees(newlon)
            mylats[icounter] = newlat_degrees
            mylons[icounter] = newlon_degrees
            icounter += 1
            f = f + fractionalincrement

    # Write the ending coordinates
    mylats.append([])
    mylons.append([])
    mylats[onelessthansegments] = ptlat2
    mylons[onelessthansegments] = ptlon2

    return mylons, mylats
