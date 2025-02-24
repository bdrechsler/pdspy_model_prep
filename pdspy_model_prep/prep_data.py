#!/usr/bin/env python
import glob
import os
import numpy as np
from casatasks import mstransform, tclean, listobs, concat, exportfits, msmd
import pdspy.interferometry as uv
from .gas_lines import line_dict
from .create_config import create_config
from .create_batch_submit import create_batch_submit

def prep_data(source_name, source_dir, line_name, chan_width, nchan, 
              vsys, robust, remove_files):
    """
    Prep the data by performing an mstranform and imaging the data
    Also perform a UV cut and eliminate the data with 0 weight
    Output the data in the correct hdf5 file format
    
    Args:
        source_name (str): The name of the source being modeled
        source_dir (str): The path to the root directory of the model directory
        line_name (str): The name of the line being modeled
        chan_width (str): The velocity with of the channels (need to include units in string) 
        nchan (int): The number of velocity channels in the image cube
        vsys (float): The system velocity of the source (in km/s)
        robust (float): The robust value used in the imaging
        remove_files (bool): If True, will delete the extra data files to save space
    """
    line_vis_list = glob.glob(source_dir + "data/*spectral_line.ms") # list of spectral line ms files
    line = line_dict[line_name]
    rfreq = line['rest_freq']

    """
    #### Regrid the data #####
    """
    regrid_vis = []

    for line_vis in line_vis_list:
        # find spectral window corresponding to rfreq
        # convert rest frequency to Hz
        rfreq_Hz = float(rfreq[:-3]) * 1e9
        
        # use msmd to get spw info
        msmd.open(line_vis)
        nspw = msmd.nspw()
        # create an array of the central freq for each spw
        spw_freqs = np.array([msmd.meanfreq(spw) for spw in range(nspw)])

        # find spw with central freq closest to rfreq 
        spw = str(np.argmin(np.abs(spw_freqs - rfreq_Hz)))
        
        outfile = line_vis.replace('spectral_line.ms', line_name+'_lsrk.ms') # name of output ms file
        print("Creating " + outfile)
        os.system("rm -rf " + outfile)

        svel = str(vsys - 7.0) + "km/s"

        mstransform(vis=line_vis, regridms=True, mode='velocity', outputvis=outfile,
                    outframe='LSRK', veltype='radio', restfreq=rfreq, datacolumn='data',
                    width=chan_width, nchan=nchan, start=svel, combinespws=False, keepflags=False,
                    timeaverage=True, timebin='30.25s', spw=spw)

        regrid_vis.append(outfile)

    """
    ##### Clean and Image the Data #######
    """
    data_dir = source_dir + 'data/'
    imagename = data_dir + source_name + '_' + line_name + '_t2000klam'
    for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.workingdirectory']:
        os.system('rm -rf ' + imagename + ext)

    uvrange='>50klambda'

    tclean(vis=regrid_vis, spw='', imagename=imagename, specmode='cube', imsize=512,
        deconvolver='hogbom', start=svel, width=chan_width, nchan=nchan, outframe='LSRK',
        veltype='radio', restfreq=rfreq, cell='0.025arcsec',
        gain=0.1, niter=20000, weighting='briggs', robust=robust, threshold='1.0mJy',
        usemask='auto-multithresh', sidelobethreshold=2.0, noisethreshold=4.0,
        lownoisethreshold=1.0, interactive=False, restoringbeam='common',
        uvtaper=['2000klambda'], uvrange=uvrange, parallel=False)

    # export the images
    myimages = glob.glob(data_dir + "*.image")
    for image in myimages:
        exportfits(imagename=image, fitsimage=image+'.fits', overwrite=True)

    for ext in ['*.pb', '*.psf', '*.model', '*.residual', '*.mask', '*.image', '*.pbcor']:
        os.system("rm -rf " + data_dir + ext)

    """
    ###### Write Data to an HDF5 file #######
    """
    print("Writting hdf5")
    # concatinate regridded visibility files
    os.system("rm -rf " + data_dir + "*concat.ms")
    concat_file = data_dir + source_name + '_' + line_name + '_' + 'concat.ms'
    concat(vis=regrid_vis, concatvis=concat_file)

    # load concated ms file into pdspy
    data = uv.readms(filename=concat_file, datacolumn='data')

    # find indicies where weights are non-zero
    good, = np.where((data.weights[:,0] > 0) & (data.uvdist > 50000))

    # find data values at these indicies
    new_u = data.u[good]
    new_v = data.v[good]
    new_real = data.real[good,:]
    new_imag = data.imag[good,:]
    new_weights = data.weights[good,:]
    # create hdf5 file
    os.system("rm -rf {}/*.hdf5".format(data_dir))
    output_file = data_dir + source_name + '_' + line_name + '_50klam.hdf5'
    new_data = uv.Visibilities(new_u, new_v, data.freq, new_real, new_imag, new_weights)
    new_data.write(output_file)
    if remove_files:
        os.system('rm -rf {}/*.ms'.format(data_dir))
        os.system('rm -rf {}/*.listobs'.format(data_dir))
        os.system('rm -rf {}/*.sumwt'.format(data_dir))
