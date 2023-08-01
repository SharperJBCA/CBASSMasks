# Example of generating region masks 
# Can define masks by thresholds, and ellipses.
# 
# SH - 28/07/2023 : Version 1. 

import numpy as np
from matplotlib import pyplot
import sys
import os
import healpy as hp 
from tools import healpix_functions
import matplotlib.colors as mcolors
from Regions import Regions, Region
from DefineRegions import threshold_mask, remove_pixels_with_masked_neighbour,query_arc_pixels, remove_overlap_pixels,smooth_mask,query_elliptical_pixels, latitude_cut, bad_data_mask, todegrees,planck_source_catalogue,subtract_elliptical_pixels

    
def create_masks(): 
    cbass_map = hp.read_map('/scratch/nas_cbassarc/cbass_data/Reductions/v34m3_mcal1/NIGHTMERID20/AWR1/calibrated_map/AWR1_xND12_xAS14_1024_NM20S3M1_G_Offmap.fits')
        #'/path/to/cbass/map/AWR1_xND12_xAS14_1024_NM20S3M1_G_Offmap.fits')
    nside = 1024
    ud_cbass_map = hp.ud_grade(cbass_map,nside)
    print(np.min(ud_cbass_map),np.max(ud_cbass_map))
    regions = Regions()

#(query_elliptical_pixels,{'nside':nside, 'ra':335.0,'dec':17.5,'major':67,'minor':67,'pa':0,'res':12})
    regions.add_region(Region(name='North Polar Spur A', cmap='Purples',nside=nside,
                                processes = [(query_arc_pixels,{'nside':nside, 'ra':329,'dec':17.5,'major':58,'minor':58,'width':30,'pa':0,'res':6,'theta_start':0,'theta_end':np.radians(50)}),                                             
                                            (subtract_elliptical_pixels,{'nside':nside, 'ra':0,'dec':0,'major':13,'minor':60,'pa':0,'res':18}), # Inner Plane
                                            (subtract_elliptical_pixels,{'nside':nside, 'ra':6,'dec':24,'major':5,'minor':5,'pa':0,'res':18}), # Sh2-27
                                            (latitude_cut,{'latitude':0,'model':'<','nside':nside}),
                                            (bad_data_mask,{'map':hp.ud_grade(cbass_map,nside)})]))
    
    regions.add_region(Region(name='North Polar Spur', cmap='Purples',nside=nside,
                                processes = [(query_arc_pixels,{'nside':nside, 'ra':329,'dec':17.5,'major':58,'minor':58,'width':30,'pa':0,'res':6,'theta_start':0,'theta_end':np.radians(160)}),                                             
                                            (subtract_elliptical_pixels,{'nside':nside, 'ra':0,'dec':0,'major':13,'minor':60,'pa':0,'res':18}), # Inner Plane
                                            (subtract_elliptical_pixels,{'nside':nside, 'ra':6,'dec':24,'major':5,'minor':5,'pa':0,'res':18}), # Sh2-27
                                            (latitude_cut,{'latitude':0,'model':'<','nside':nside}),
                                            (bad_data_mask,{'map':hp.ud_grade(cbass_map,nside)})]))
    regions.add_region(Region(name='North Polar Spur B', cmap='Purples',nside=nside,
                                processes = [(query_arc_pixels,{'nside':nside, 'ra':329,'dec':17.5,'major':58,'minor':58,'width':30,'pa':0,'res':6,'theta_start':np.radians(50),'theta_end':np.radians(120)}),                                             
                                            (subtract_elliptical_pixels,{'nside':nside, 'ra':0,'dec':0,'major':13,'minor':60,'pa':0,'res':18}), # Inner Plane
                                            (subtract_elliptical_pixels,{'nside':nside, 'ra':6,'dec':24,'major':5,'minor':5,'pa':0,'res':18}), # Sh2-27
                                            (latitude_cut,{'latitude':0,'model':'<','nside':nside}),
                                            (bad_data_mask,{'map':hp.ud_grade(cbass_map,nside)})]))

    regions.add_region(Region(name='CBASS 95pc', cmap='Reds',nside=nside,
                              processes=[(threshold_mask,{'map':ud_cbass_map,'threshold_pc':95}),
                                         (remove_pixels_with_masked_neighbour,{'nside':nside}),
                                         (smooth_mask,{'nside':nside,'fwhm':1}),
                                         (bad_data_mask,{'map':ud_cbass_map})]))
    regions.add_region(Region(name='CBASS 90pc', cmap='Oranges',nside=nside,
                              processes=[(threshold_mask,{'map':ud_cbass_map,'threshold_pc':90}),
                                         (remove_pixels_with_masked_neighbour,{'nside':nside}),
                                         (smooth_mask,{'nside':nside,'fwhm':1}),                                        
                                         (bad_data_mask,{'map':ud_cbass_map})]))
    regions.add_region(Region(name='CBASS 85pc', cmap='Blues',nside=nside,
                              processes=[(threshold_mask,{'map':ud_cbass_map,'threshold_pc':85}),
                                         (remove_pixels_with_masked_neighbour,{'nside':nside}),
                                         (smooth_mask,{'nside':nside,'fwhm':1}),
                                         (bad_data_mask,{'map':ud_cbass_map})]))
    regions.add_region(Region(name='CBASS 80pc', cmap='Blues',nside=nside,
                              processes=[(threshold_mask,{'map':ud_cbass_map,'threshold_pc':80}),
                                         (remove_pixels_with_masked_neighbour,{'nside':nside}),
                                         (smooth_mask,{'nside':nside,'fwhm':1}),
                                         (bad_data_mask,{'map':ud_cbass_map})]))

    regions.plot_regions(background=cbass_map, figname='percent_cuts.png', cmap=pyplot.cm.viridis)

    regions.add_region(Region(name='LFI30 PCCS', cmap='Blues',nside=nside,
                              processes=[(planck_source_catalogue,dict(nside=128, catalogue='../other_masks/source_catalogues/lfi30_pccs_planck2016.fits', hdu_index=1, flux_threshold=1,min_latitude=10)),
                                         (bad_data_mask,{'map':hp.ud_grade(cbass_map,nside)})]))

    regions.add_region(Region(name='Orion-Eridanus superbubble', cmap='Blues',nside=nside,
                              processes = [(query_elliptical_pixels,{'nside':nside, 'ra':195,'dec':-30,'major':22.5,'minor':10,'pa':-45,'res':18})]))
    regions.add_region(Region(name='Fan Region', cmap='Blues',nside=nside,
                                processes = [(query_elliptical_pixels,{'nside':nside, 'ra':135,'dec':0,'major':30,'minor':40,'pa':0,'res':18}),
                                            (threshold_mask,{'map':hp.ud_grade(cbass_map,nside),'threshold':-0.005}),
                                            (remove_pixels_with_masked_neighbour,{'nside':nside}),
                                            (remove_pixels_with_masked_neighbour,{'nside':nside}),
                                            (smooth_mask,{'nside':nside,'fwhm':1})]))
    regions.add_region(Region(name='Inner Plane', cmap='Blues',nside=nside,
                                processes = [(query_elliptical_pixels,{'nside':nside, 'ra':0,'dec':0,'major':13,'minor':60,'pa':0,'res':18})]))
    regions.add_region(Region(name='Sh2-27', cmap='Blues',nside=nside,
                                processes = [(query_elliptical_pixels,{'nside':nside, 'ra':6,'dec':24,'major':5,'minor':5,'pa':0,'res':18})]))    
    regions.add_region(Region(name='Loop II', cmap='Reds',nside=nside,
                                processes = [(query_arc_pixels,{'nside':nside, 'ra':100,'dec':-32,'major':45.5,'minor':45.5,'width':20,'pa':0,'res':6,'theta_start':np.radians(200),'theta_end':np.radians(330)}),
                                            (bad_data_mask,{'map':hp.ud_grade(cbass_map,nside)})]))
    regions.add_region(Region(name='Loop III', cmap='Reds',nside=nside,
                                processes = [(query_arc_pixels,{'nside':nside, 'ra':124,'dec':15.5,'major':32.5,'minor':32.5,'width':20,'pa':0,'res':6,'theta_start':np.radians(20),'theta_end':np.radians(170)}),
                                             (latitude_cut,{'latitude':0,'model':'<','nside':nside}),
                                            (bad_data_mask,{'map':hp.ud_grade(cbass_map,nside)})]))
    regions.add_region(Region(name='Loop IIIs', cmap='Reds',nside=nside,
                                processes = [(query_arc_pixels,{'nside':nside, 'ra':106,'dec':-22,'major':25,'minor':25,'width':15,'pa':0,'res':6,'theta_start':np.radians(200),'theta_end':np.radians(330)}),
                                            (bad_data_mask,{'map':hp.ud_grade(cbass_map,nside)})]))

    regions.add_region(Region(name='R3', cmap='Greens',nside=nside,
                                processes = [(query_elliptical_pixels,{'nside':nside, 'ra':110,'dec':60,'major':15,'minor':15,'pa':0,'res':18})]))
    regions.add_region(Region(name='Minimum', cmap='gray',nside=nside,
                                processes = [(query_elliptical_pixels,{'nside':nside, 'ra':210,'dec':45,'major':30,'minor':30,'pa':0,'res':18})]))


    regions.plot_regions(background=cbass_map, figname='test.png', cmap=pyplot.cm.viridis)

    regions.write_regions('regions.h5',overwrite=True)

def plot_masks():
    regions = Regions()
    regions.load_regions('regions.h5')

    for nside in [1024,512,256, 128]:
        regions.create_mask(f'masks/CBASS_95pc_G_{nside:04d}.fits',include=['CBASS 95pc'],nside=nside)
        regions.create_mask(f'masks/CBASS_90pc_G_{nside:04d}.fits',include=['CBASS 90pc'],nside=nside)
        regions.create_mask(f'masks/CBASS_85pc_G_{nside:04d}.fits',include=['CBASS 85pc'],nside=nside)
        regions.create_mask(f'masks/CBASS_80pc_G_{nside:04d}.fits',include=['CBASS 80pc'],nside=nside)
        regions.create_mask(f'masks/Fan_Region_G_{nside:04d}.fits',include=['Fan Region'],nside=nside)
        regions.create_mask(f'masks/NorthPolarSpur_G_{nside:04d}.fits',include=['North Polar Spur'],exclude=['LFI30 PCCS'],nside=nside)
        regions.create_mask(f'masks/NorthPolarSpurA_G_{nside:04d}.fits',include=['North Polar Spur A'],exclude=['LFI30 PCCS'],nside=nside)
        regions.create_mask(f'masks/NorthPolarSpurB_G_{nside:04d}.fits',include=['North Polar Spur B'],exclude=['LFI30 PCCS'],nside=nside)


if __name__ == "__main__":
    create_masks()
    #plot_masks()
