import numpy as np
from matplotlib import pyplot
import sys
import os 
from dataclasses import dataclass, field
from tqdm import tqdm
import h5py 
from tools import healpix_functions
import matplotlib.colors as mcolors
from astropy.io import fits 
from tools.healpix_functions import todegrees, tothetaphi
import healpy as hp
import matplotlib.patheffects as pe

@dataclass
class Region: 
    name : str = 'none'
    cmap : str = 'viridis'
    nside : int = 128
    processes : dict = field(default_factory=dict)
    pixels : np.ndarray = field(default_factory=lambda:np.zeros(1))

    def __post_init__(self):
        print(self.name)
        self.pixels = np.arange(12*self.nside**2,dtype=int)
        for process, parameters in self.processes:
            self.pixels = process(self.pixels,**parameters)
        if 'pc' in self.name:
            print(self.pixels)
    def __items__(self):
        return self.processes.items()
    
    def __keys__(self):
        return self.processes.keys()
    
    def __values__(self):
        return self.processes.values()

class Regions: 
    def __init__(self):
        self.regions = []

    def list_masks(self):
        for region in self.regions:
            print(region.name)
        return self.mask_names
    
    @property
    def mask_names(self):
        return [region.name for region in self.regions]
    
    def get_regions(self, names):
        return [region for region in self.regions if region.name in names]

    def add_region(self, region : Region):
        self.regions.append(region)

    def load_regions(self, filename : str):
        h = h5py.File(filename,'r')
        for region_name, region_info in h.items():
            region = Region()
            region.name = region_name
            region.pixels = np.array(region_info['pixels'][...])
            region.nside = region_info.attrs['nside']
            region.cmap = region_info.attrs['cmap']
            self.add_region(region)
        h.close()

    def write_regions(self, filename : str, overwrite : bool=False):
        if os.path.exists(filename) and overwrite:
            os.remove(filename)
        h = h5py.File(filename,'w')
        for region in self.regions:
            region_info = h.create_group(region.name)
            region_info.create_dataset('pixels',data=region.pixels)
            region_info.attrs['nside'] = region.nside
            region_info.attrs['cmap'] = region.cmap
        h.close()

    def create_mask(self, mask_filename, include=[], exclude=[]):
        if not os.path.exists(os.path.dirname(mask_filename)):
            os.makedirs(os.path.dirname(mask_filename))

        # First get the include pixels
        include_pixels = []
        for region in self.regions:
            if region.name in include:
                include_pixels.append(region.pixels)
        include_pixels = np.concatenate(include_pixels)

        # Then get the exclude pixels
        exclude_pixels = []
        for region in self.regions:
            if region.name in exclude:
                exclude_pixels.append(region.pixels)
        if len(exclude_pixels) > 0:
            exclude_pixels = np.concatenate(exclude_pixels)
            # Remove exclude pixels from include pixels
            nside = self.regions[0].nside
            m = np.zeros(12*nside**2)
            m[include_pixels] = 1
            m[exclude_pixels] += 2
            hp.mollview(m)
            pyplot.savefig('test.png')
            pyplot.close()
            mask_pixels = np.where(m == 1)[0]
        else:
            mask_pixels = include_pixels

        # Create the mask
        mask = np.zeros(12*self.regions[0].nside**2)
        mask[mask_pixels] = 1

        # Write the mask
        hp.write_map(mask_filename,mask,overwrite=True)

        tmpregion = Region(name='',pixels=mask_pixels,nside=self.regions[0].nside,cmap=self.get_regions(include)[0].cmap)
        self.plot_regions([tmpregion],background=mask,figname=mask_filename.replace('.fits','.png'))

    def plot_regions(self, regions=None, background=None, figname='regions.png', cmap=pyplot.cm.viridis):

        # Convert colormap to 256 RGBA values
        colors = cmap(np.linspace(0., 1., cmap.N))

        # Make the alpha vary from 0 to 1
        colors[:,-1] = np.linspace(0, 1, cmap.N)

        # Create a new colormap
        cmap_alpha = mcolors.ListedColormap(colors)

        mollview = healpix_functions.Mollview()

        fig = pyplot.figure(figsize=(10,10))
        if isinstance(background,type(None)):
            background = np.zeros(12*self.regions[0].nside**2)
            cmap = pyplot.cm.get_cmap('Greys')
        mollview(background,vmin=-0.05,vmax=0.05,cmap=cmap)

        if regions is None:
            regions = self.regions
        for region in regions:
            m = np.zeros(12*region.nside**2)
            m[region.pixels] = 1
            m[m==0]=np.nan
            mollview.contourf(m,cmap=pyplot.get_cmap(region.cmap))
            mollview.text(*todegrees(*hp.pix2ang(region.nside,region.pixels[0])),region.name,fontsize=10,color='w',path_effects=[pe.withStroke(linewidth=1, foreground="k")])
        #hp.graticule()
        pyplot.savefig(figname,dpi=300)
