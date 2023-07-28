import numpy as np
from matplotlib import pyplot
import sys
import os 
from dataclasses import dataclass, field
from tqdm import tqdm
import h5py 

@dataclass
class Region:

    name : str = '' 
    pixels : np.ndarray = field(default_factory=np.array)
    nside : int = 0 
    cmap : str = 'viridis'


class Regions: 
    def __init__(self):
        self.regions = []

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

    def write_regions(self, filename : str):
        h = h5py.File(filename,'w')
        for region in self.regions:
            region_info = h.create_group(region.name)
            region_info.create_dataset('pixels',data=region.pixels)
            region_info.attrs['nside'] = region.nside
            region_info.attrs['cmap'] = region.cmap
        h.close()