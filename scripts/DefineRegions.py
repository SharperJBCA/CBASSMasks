import numpy as np
from matplotlib import pyplot
import sys
import os
import healpy as hp 
from tools import healpix_functions
import matplotlib.colors as mcolors
from Regions import Regions, Region

def query_elliptical_pixels(nside, ra, dec, major, minor=None, pa=0, radians_input = False,theta_start=0,theta_end=2*np.pi,res=12):
    """Query pixels within an elliptical region.

    Parameters
    ----------
    nside : int
        Healpix nside.
    ra : float
        Right ascension of the center of the ellipse (degrees).
    dec : float
        Declination of the center of the ellipse (degrees).
    major : float
        Major axis of the ellipse (degrees).
    minor : float
        Minor axis of the ellipse (degrees).
    pa : float
        Position angle of the ellipse (degrees).

    Returns
    -------
    pixels : array-like
        Healpix pixels within the ellipse.
    """
    # Convert to radians
    if minor is None:
        minor = major

    if not radians_input:
        ra = np.radians(ra)
        dec = np.radians(dec)
        major = np.radians(major)
        minor = np.radians(minor)
        pa = np.radians(pa)

    # Define the ellipse
    theta = np.linspace(theta_start, theta_end, res, endpoint=False)
    x = np.cos(theta) * minor
    y = np.pi/2. - np.sin(theta) * major

    # Now rotate the ellipse to the ra, dec, pa location 
    rot = hp.rotator.Rotator(rot=[np.degrees(ra),np.degrees(dec),np.degrees(pa)], inv=True) 
    theta, phi = rot(y, x)

    # Find all pixels bound by ellipse 
    pixels = hp.query_polygon(nside, hp.ang2vec(theta[::-1], phi[::-1]))
    
    return pixels

def threshold_mask(pixels, map, threshold):
    """Threshold a map and return the pixels above the threshold.

    Parameters
    ----------
    pixels : array-like
        Healpix pixels to threshold.
    map : array-like
        Healpix map to threshold.
    threshold : float
        Threshold value.

    Returns
    -------
    pixels : array-like
        Healpix pixels above threshold.
    """
    return pixels[np.where(map[pixels] > threshold)[0]]

def remove_pixels_with_masked_neighbour(nside,pixels):
    """Remove pixels that have a masked neighbour.

    Parameters
    ----------
    pixels : array-like
        Healpix pixels to threshold.

    Returns
    -------
    pixels : array-like
        Healpix pixels above threshold.
    """
    # Get the neighbours of the pixels
    m = np.zeros(12*nside**2)
    m[pixels] = 1

    neighbours = hp.get_all_neighbours(nside, pixels)

    # Remove the pixels that have a masked neighbour
    print(np.sum(m[neighbours], axis=0))
    return pixels[np.where(np.sum(m[neighbours], axis=0) > 6)[0]]

def remove_overlap_pixels(nside, pixels_main, pixels_cut):

    m = np.zeros(12*nside**2)
    m[pixels_main] = 1
    m[pixels_cut] += 2

    return np.where(m == 1)[0]

def latitude_cut(nside, pixels, latitude, model='>'):
    """Cut pixels based on latitude.

    Parameters
    ----------
    pixels : array-like
        Healpix pixels to threshold.
    latitude : float
        Latitude to cut at (degrees).
    model : str
        Cut model. Can be '>' or '<'.

    Returns
    -------
    pixels : array-like
        Healpix pixels above threshold.
    """
    theta, phi = hp.pix2ang(nside, pixels)
    theta = np.pi/2. - theta
    latitude = np.radians(latitude)
    if model == '<':
        return pixels[np.where(theta > latitude)[0]]
    elif model == '>':
        return pixels[np.where(theta < latitude)[0]]
    else:
        raise ValueError('Invalid model.')
    
def bad_data_mask(m,pixels):
    """Mask bad data.

    Parameters
    ----------
    m : array-like
        Healpix map.
    pixels : array-like
        Healpix pixels to threshold.

    Returns
    -------
    pixels : array-like
        Healpix pixels above threshold.
    """
    return pixels[(m[pixels] != hp.UNSEEN)]

def todegrees(theta,phi):
    return np.degrees(phi), np.degrees(np.pi/2-theta)

if __name__ == "__main__":
    cbass_map = hp.read_map('/scratch/nas_cbassarc/cbass_data/Reductions/v34m3_mcal1/NIGHTMERID20/AWR1/calibrated_map/AWR1_xND12_xAS14_1024_NM20S3M1_G_Offmap.fits')
    nside = 128
    Regions = Regions()
    # Orion-Eridanus superbubble
    ra = 195
    dec = -30
    major = 22.5
    minor = 10 
    pa = -45
    pixels = query_elliptical_pixels(nside, ra, dec, major, minor=minor, pa=pa, radians_input = False)
    Regions.add_region(Region('Orion-Eridanus superbubble',pixels,nside,'Blues'))

    # Fan Region
    ra = 135
    dec = 0
    major = 30
    minor = 40 
    pa = 0 
    pixels = query_elliptical_pixels(nside, ra, dec, major, minor=minor, pa=pa, radians_input = False)
    pixels = threshold_mask(pixels,hp.ud_grade(cbass_map,nside),-0.005)
    pixels = remove_pixels_with_masked_neighbour(nside,pixels)
    pixels = remove_pixels_with_masked_neighbour(nside,pixels)
    print(pixels.size)
    Regions.add_region(Region('Fan Region',pixels,nside,'Blues'))

    # Inner Plane
    ra = 0
    dec = 0
    major = 13
    minor = 60
    pa = 0
    pixels_plane = query_elliptical_pixels(nside, ra, dec, major, minor=minor, pa=pa, radians_input = False)
    Regions.add_region(Region('Inner Plane',pixels_plane,nside,'Blues'))

    # Hii Region
    ra = 6
    dec = 24
    major = 5
    minor = 5
    pa = 0
    pixels_hii = query_elliptical_pixels(nside, ra, dec, major, minor=minor, pa=pa, radians_input = False)
    Regions.add_region(Region('Sh2-27',pixels_hii,nside,'Blues'))

    # North Polar Spur
    ra = 335.0
    dec = 17.5
    major = 67
    minor = major
    pa = 0
    pixels_outer = query_elliptical_pixels(nside, ra, dec, major, minor=minor, pa=pa, radians_input = False)
    # Inner North Polar Spur
    ra = 329.0
    dec = 17.5
    major = 40
    minor = major
    pa = 0
    pixels_inner = query_elliptical_pixels(nside, ra, dec, major, minor=minor, pa=pa, radians_input = False)

    pixels = remove_overlap_pixels(nside, pixels_outer,pixels_inner) 
    pixels = remove_overlap_pixels(nside, pixels,pixels_plane) 
    pixels = remove_overlap_pixels(nside, pixels,pixels_hii) 
    pixels = latitude_cut(nside,pixels,0, model='<')
    pixels = bad_data_mask(hp.ud_grade(cbass_map,nside),pixels) 
    Regions.add_region(Region('North Polar Spur',pixels,nside,'Reds'))

    # Loop II
    ra = 100
    dec = -32.5
    major = 55.5
    minor = major
    pa = 0
    theta_start = np.pi
    theta_end = 2*np.pi
    pixels_outer = query_elliptical_pixels(nside, ra, dec, major, minor=minor, pa=pa, radians_input = False,theta_start=theta_start,theta_end=theta_end)
    # Inner North Polar Spur
    ra = 100
    dec = -32.5
    major = 35.5
    minor = major
    pa = 0
    theta_start = np.pi
    theta_end = 2*np.pi
    pixels_inner = query_elliptical_pixels(nside, ra, dec, major, minor=minor, pa=pa, radians_input = False,theta_start=theta_start,theta_end=theta_end)

    pixels = remove_overlap_pixels(nside, pixels_outer,pixels_inner) 
    pixels = bad_data_mask(hp.ud_grade(cbass_map,nside),pixels) 

    Regions.add_region(Region('Loop II',pixels,nside,'Reds'))

    # Loop III
    ra = 124
    dec = 15.5
    major = 37.5
    minor = major
    pa = 0
    theta_start = np.radians(20)
    theta_end = np.pi
    pixels_outer = query_elliptical_pixels(nside, ra, dec, major, minor=minor, pa=pa, radians_input = False,theta_start=theta_start,theta_end=theta_end)
    # Inner 
    ra = 124
    dec = 15.5
    major = 28.5
    minor = major
    pa = 0
    theta_start = np.radians(20)
    theta_end = np.pi
    pixels_inner = query_elliptical_pixels(nside, ra, dec, major, minor=minor, pa=pa, radians_input = False,theta_start=theta_start,theta_end=theta_end)

    pixels = remove_overlap_pixels(nside, pixels_outer,pixels_inner) 
    pixels = bad_data_mask(hp.ud_grade(cbass_map,nside),pixels) 

    Regions.add_region(Region('Loop III',pixels,nside,'Reds'))

    # R3
    ra = 110
    dec = 60
    major = 15
    minor = major
    pa = 0
    theta_start = 0
    theta_end = 2*np.pi
    pixels = query_elliptical_pixels(nside, ra, dec, major, minor=minor, pa=pa, radians_input = False,theta_start=theta_start,theta_end=theta_end,res=4)
    Regions.add_region(Region('R3',pixels,nside,'Greens'))

    # Minimum
    ra = 210
    dec = 45
    major = 30
    minor = major
    pa = 0
    theta_start = 0
    theta_end = 2*np.pi
    pixels = query_elliptical_pixels(nside, ra, dec, major, minor=minor, pa=pa, radians_input = False,theta_start=theta_start,theta_end=theta_end,res=4)
    Regions.add_region(Region('Minimum',pixels,nside,'gray'))





    # Let's use the viridis colormap as an example
    cmap = pyplot.cm.viridis

    # Convert colormap to 256 RGBA values
    colors = cmap(np.linspace(0., 1., cmap.N))

    # Make the alpha vary from 0 to 1
    colors[:,-1] = np.linspace(0, 1, cmap.N)

    # Create a new colormap
    cmap_alpha = mcolors.ListedColormap(colors)

    mollview = healpix_functions.Mollview()

    fig = pyplot.figure(figsize=(10,10))
    mollview(cbass_map,vmin=-0.05,vmax=0.05,cmap=pyplot.cm.viridis)

    for region in Regions.regions:
        m = np.zeros(12*nside**2)
        m[region.pixels] = 1
        m[m==0]=np.nan
        mollview.contourf(m,cmap=pyplot.get_cmap(region.cmap))
        mollview.text(*todegrees(*hp.pix2ang(nside,region.pixels[0])),region.name,fontsize=10,color='w')
    #hp.graticule()
    pyplot.savefig('test.png')

    Regions.write_regions('regions.h5')