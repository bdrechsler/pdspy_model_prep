import os
import glob
from .prep_data import prep_data
from .create_config import create_config
from .create_batch_submit import create_batch_submit

class Model:

    def __init__(self, source_name, dpc, vsys, line_name='C18O', chan_width='0.334km/s',
                 nchan=42, x0=[-0.5, 0.5], y0=[-0.5, 0.5], source_dir = './', svel=None,
                 disk_types=["truncated", "exptaper", "dartois-exptaper", "dartois-truncated"]):

        self.source_name = source_name
        self.dpc = dpc
        self.line_name = line_name
        self.chan_width = chan_width
        self.nchan = nchan
        self.vsys = vsys
        self.svel = svel
        self.x0=x0
        self.y0=y0
        self.source_dir = source_dir
        self.disk_types = disk_types

    def prep_model(self, data=True, config=True, batch_script=True, robust=2.0,
                   remove_files=False, ncpu=1):

        # if there are ms files, prep the data
        data_dir = self.source_dir + 'data/'
        ms_files = glob.glob(data_dir + "*.ms")
        if data and len(ms_files) != 0:
            
            prep_data(source_name = self.source_name, source_dir = self.source_dir, 
                      line_name = self.line_name, chan_width = self.chan_width, 
                      nchan = self.nchan, vsys = self.vsys, robust= robust,
                      remove_files = remove_files, svel = self.svel)

        if config or batch_script:
            for disk_type in self.disk_types:
                model_dir = self.source_dir + disk_type + '_' + self.line_name + '/'
                if not os.path.exists(model_dir):
                    os.system("mkdir {}".format(model_dir))

        if config:
            for disk_type in self.disk_types:
                create_config(source_name = self.source_name, source_dir=self.source_dir, 
                              line_name = self.line_name, disk_type = disk_type,
                              dpc = self.dpc, vsys = self.vsys, x0 = self.x0, y0 = self.y0)

        if batch_script:
            for disk_type in self.disk_types:
                create_batch_submit(source_name = self.source_name, line_name=self.line_name, 
                                    source_dir = self.source_dir, disk_type = disk_type,
                                    ncpu = ncpu)
