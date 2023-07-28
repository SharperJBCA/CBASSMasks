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

