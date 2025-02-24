import numpy as np
from .gas_lines import line_dict

def write_dict(d, d_name, file):
    """
    Write a dictionary to a file, meant to write a dictionary to the 
    config.py file
    Args:
        d (dict): dictionary to write
        d_name (string): name of the dictionary in the new file
        file (string): path to the file to append the dictionary to
    """
    file.write(d_name + ' = {\n')
    for key, val in d.items():
        file.write('        "{0}":{1},\n'.format(key, val).replace('array', 'np.array'))
    file.write('        }\n\n')


def write_dynesty_params(file, nwalkers=100, steps_per_iter=20, max_nsteps=2000,
                        nplot=5, nlive_init=250, nlive_batch=250, maxbatch=0,
                        dlogz=0.05, walks=25):
    """
    Write the standard dynesty parameters to the config.py file

    """
    
    file.write('nwalkers = {}\n'.format(nwalkers))
    file.write('steps_per_iter = {}\n'.format(steps_per_iter))
    file.write('max_nsteps = {}\n'.format(max_nsteps))
    file.write('nplot = {}\n'.format(nplot))
    file.write('nlive_init = {}\n'.format(nlive_init))
    file.write('nlive_batch = {}\n'.format(nlive_batch))
    file.write('maxbatch = {}\n'.format(maxbatch))
    file.write('dlogz = {}\n'.format(dlogz))
    file.write('walks = {}\n\n'.format(walks))

def create_config(source_name, source_dir, line_name, disk_type, dpc, vsys, x0, y0):
    """Create a config.py file for the fit run"""
    
    source_dir = source_dir + source_name + '/'
    # paths to the data (hdf5) and image (fits) files
    data_file = "../data/{0}_{1}_50klam.hdf5".format(source_name, line_name)
    image_file= "../data/{0}_{1}_t2000klam.image.fits".format(source_name, line_name)

    line = line_dict['line_name']
    freq = line['rest_freq']
    gas_file = line['data_file']

    visibilities = {
            "file":[data_file],
            "binsize":[8057.218995847603],
            "pixelsize":[0.01],
            "freq":[freq],
            "lam":[line_name],
            "npix":[1024],
            "x0":[0.0],
            "y0":[0.0],
            "gridsize":[1024],
            "weight":[0.25],
            # Info for the image.
            "image_file":[image_file],
            "image_pixelsize":[0.025],
            "image_npix":[512],
            # Info for the plots.
            "nrows":[7],
            "ncols":[8],
            "ind0":[1],
            "fmt":['5.2f'],
            "ticks":[np.array([-250,-200,-100,0,100,200,250])],
            "image_ticks":[np.array([-0.75,-0.5,0,0.5,0.75])],
            "nphi":[256],
            "nr":[2],
            }

    spectra = {
            "file":[],
            "bin?":[],
            "nbins":[],        
            "weight":[],
            }


    images = {
            "file":[],
            "npix":[],
            "pixelsize":[],
            "lam":[],
            "bmaj":[],
            "bmin":[],
            "bpa":[],
            "ticks":[],
            "plot_mode":[],
            }


    parameters = {
            # Stellar parameters.
            "logM_star":{"fixed":False, "value":-0.3, "limits":[-1.5,1.0]},
            # Disk parameters.
            "disk_type":{"fixed":True, "value":disk_type, "limits":[0.,0.]},
            "logM_disk":{"fixed":False, "value":-3.0, "limits":[-10.,-2.5]},
            "logR_disk":{"fixed":False, "value":2, "limits":[1.,2.6]},
            "gamma":{"fixed":False, "value":1.0, "limits":[0.0,2.0]},
            "f_M_large":{"fixed":True, "value":0.8, "limits":[0.05, 1.0]},
            # Disk temperature parameters.
            "logT0":{"fixed":False, "value":2.5, "limits":[1.,3.0]},
            "q":{"fixed":True, "value":0.25, "limits":[0.,1.]},
            "loga_turb":{"fixed":True, "value":-3.0, "limits":[-3.,-0.875]},
            # Envelope parameters.
            "envelope_type":{"fixed":True, "value":"ulrich-tapered", "limits":[0.,0.]},
            "logM_env":{"fixed":False, "value":-3., "limits":[-6., -2.]},
            "logR_env":{"fixed":False, "value":4., "limits": [2.,4.]},
            "loga_turb_env":{"fixed":True, "value":-3.0, "limits":[-3.,-0.875]},
            # Gas parameters.
            "gas_file1":{"fixed":True, "value":gas_file, "limits":[0.,0.]},
            "logabundance1":{"fixed":True, "value":np.log10(1.0e-4/560.), "limits":[-10.,-2.]},
        # Viewing parameters.
            "i":{"fixed":False, "value":90.0, "limits":[0.,180.]},
            "pa":{"fixed":False, "value":94.85, "limits":[0.0, 360.0]},
            "x0":{"fixed":False, "value":0, "limits": x0},
            "y0":{"fixed":False, "value":0, "limits": y0},
            "dpc":{"fixed":True, "value":dpc, "limits":[1.,1e6]},
            "Ak":{"fixed":True, "value":0., "limits":[0.,1.]},
            "v_sys":{"fixed":False, "value":6.0, "limits":[vsys - 3.0, vsys + 3.0]}
            }

    # add dartois specific parameters
    if "dartois" in disk_type:
            del parameters["logT0"]
            del parameters["q"]
            parameters["logTmid0"] = {"fixed":False, "value":2.5, "limits":[0.,3.0]}
            parameters["logTatm0"] = {"fixed":False, "value":2.5, "limits":[1.,3.5]}
            parameters["pltgas_mid"] = {"fixed":False, "value":0.5, "limits":[0.,1.]}
            parameters["pltgas_atm"] = {"fixed":False, "value":"pltgas_mid", "limits":[0.,1.]}
            parameters["zq0"] = {"fixed":False, "value":0.1, "limits":[0.01,0.5]}
            parameters["plz"] = {"fixed":True, "value":1.3, "limits":[1., 1.5]}
    
    model_dir = "{0}{1}/".format(source_dir, disk_type)
    with open(model_dir + 'config.py', 'w') as config:
        # write import statements
        config.write('import numpy as np\n')
        config.write('import os\n\n')
        
        # write visibility dictionary
        write_dict(visibilities, 'visibilities', config)

        # write spectra dictionary
        write_dict(spectra, 'spectra', config)

        # write images dictionary
        write_dict(images, 'images', config)

        # write dynesty params
        write_dynesty_params(config)

        # write model params
        write_dict(parameters, 'parameters', config)
