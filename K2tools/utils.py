#!/usr/bin/evn python

import numpy as np
import matplotlib.pyplot as pl
import os

import pandas as pd
from astropy.io import fits

from matplotlib.colors import LogNorm
from scipy.ndimage import measurements
from astropy.visualization import ZScaleInterval
from astropy import table
from astropy.io import fits
from lightkurve import KeplerTargetPixelFile

#---------------------------PHOTOMETRY---------------------------#

cmap='viridis'

def make_aperture_outline(frame, no_combined_images=1, threshold=0.5):
    ## this is a little module that defines so called outlines to be used for plotting apertures

    thres_val = no_combined_images * threshold
    mapimg = (frame > thres_val)
    ver_seg = np.where(mapimg[:,1:] != mapimg[:,:-1])
    hor_seg = np.where(mapimg[1:,:] != mapimg[:-1,:])

    l = []
    for p in zip(*hor_seg):
        l.append((p[1], p[0]+1))
        l.append((p[1]+1, p[0]+1))
        l.append((np.nan,np.nan))

    # and the same for vertical segments
    for p in zip(*ver_seg):
        l.append((p[1]+1, p[0]))
        l.append((p[1]+1, p[0]+1))
        l.append((np.nan, np.nan))


    segments = np.array(l)

    x0 = -0.5
    x1 = frame.shape[1]+x0
    y0 = -0.5
    y1 = frame.shape[0]+y0

    #   now we need to know something about the image which is shown
    #   at this point let's assume it has extents (x0, y0)..(x1,y1) on the axis
    #   drawn with origin='lower'
    # with this information we can rescale our points
    segments[:,0] = x0 + (x1-x0) * segments[:,0] / mapimg.shape[1]
    segments[:,1] = y0 + (y1-y0) * segments[:,1] / mapimg.shape[0]

    return segments

def find_aperture(fluxes,starname='',kepmag='na',cutoff_limit=1.,showfig=None):
    #
    # This definition reads a 2D array of fluxes (over time) and creates an aperture mask which can later be used to select those pixels for inclusion in light curve
    #

    # first sum all the flux over the different times, this assumes limited movement throughout the time series
    flux = np.nansum(fluxes,axis=0)

    # define which cutoff flux to use for including pixel in mask
    cutoff = cutoff_limit*np.median(flux) # perhaps a more elaborate way to define this could be found in the future but this seems to work pretty well.

    # define the aperture based on cutoff and make it into array of 1 and 0
    aperture =  np.array([flux > cutoff]) #scipy.zeros((np.shape(flux)[0],np.shape(flux)[1]), int)
    aperture = np.array(1*aperture)
    #print aperture
    outline_all = make_aperture_outline(aperture[0]) # an outline (ONLY for figure) of what we are including if we would make no breakups

    # this cool little trick allows us to measure distinct blocks of apertures, and only select the biggest one
    lw, num = measurements.label(aperture) # this numbers the different apertures distinctly
    area = measurements.sum(aperture, lw, index=np.arange(lw.max() + 1)) # this measures the size of the apertures
    aperture = area[lw].astype(int) # this replaces the 1s by the size of the aperture
    aperture = (aperture >= np.max(aperture))*1 # remake into 0s and 1s but only keep the largest aperture

    outline = make_aperture_outline(aperture[0]) # a new outline (ONLY for figure)

    if showfig: # make aperture figure
        pl.figure('Aperture_' + str(starname))
        pl.imshow(flux,norm=LogNorm(),interpolation="none",cmap=cmap)
        pl.plot(outline_all[:, 0], outline_all[:, 1],color='green', zorder=10, lw=2.5)
        pl.plot(outline[:, 0], outline[:, 1],color='red', zorder=10, lw=2.5)#,label=str(kepmag))
        pl.colorbar(orientation='vertical')
        pl.xlabel('X',fontsize=15)
        pl.ylabel('Y',fontsize=15)
        #pl.legend()
        pl.tight_layout()
        pl.show()
    return np.array(aperture[0],dtype=bool)

# import warnings
# def get_centroids(flux, column, row, aperture_mask):
#     """Returns centroids based on sample moments.
#
#     Parameters
#     ----------
#     aperture_mask : array-like
#         A boolean array describing the aperture such that `False` means
#         that the pixel will be masked out.
#
#     Returns
#     -------
#     col_centr, row_centr : tuple
#         Arrays containing centroids for column and row at each cadence
#     """
#     yy, xx = np.indices(flux.shape[1:]) + 0.5
#     yy = row + yy
#     xx = column + xx
#     total_flux = np.nansum(flux[:, aperture_mask], axis=1)
#     with warnings.catch_warnings():
#         # RuntimeWarnings may occur below if total_flux contains zeros
#         warnings.simplefilter("ignore", RuntimeWarning)
#         col_centr = np.nansum(xx * aperture_mask * flux, axis=(1, 2)) / total_flux
#         row_centr = np.nansum(yy * aperture_mask * flux, axis=(1, 2)) / total_flux
#
#     return col_centr, row_centr

from photutils import centroid_com
from sklearn.ensemble import IsolationForest

def get_centroids(fluxes,  centroid_shift=None, method='abs_distance',
                    check_outliers=False, showfig=False):
    '''
    fluxes: array, 3D fluxes
    method: str, (1) absolute distance (2) sklearn.ensemble
    centroid_shift: int, maximum shift in pixels; if larger then mask is flagged
    check_outliers: bool,
    '''
    xcen, ycen = [],[]
    mask = np.zeros(len(fluxes))

    if check_outliers:
        stacked = np.nansum(fluxes, axis=0)
        Xc, Yc = centroid_com(stacked)

    for idx,img in enumerate(fluxes):
        x,y = centroid_com(img)
        if check_outliers:
            if centroid_shift is not None:
                #shift = int(round(centroid_shift))
                shift = float(centroid_shift)
                if abs(Xc-x) > shift or abs(Yc-y) > shift:
                    mask[idx] = 1
            else:
                raise ValueError('`centroid_shift` not defined')
        else:
            pass
        xcen.append(x)
        ycen.append(y)

    #convert 0,1 into bool
    mask = np.array(mask, dtype=bool)

    if method == 'ensemble':
        X = np.c_[xcen,ycen]
        mask = IsolationForest(max_features=2, contamination=0.01).fit(X).predict(X)
        #output is 1 or -1, convert to bool
        mask = np.array([True if bool(i == 1.0) else False for i in mask])

    xcen, ycen = np.array(xcen),np.array(ycen)
    if showfig:
        if method == 'ensemble':
            pl.scatter(xcen, ycen, c=mask, cmap='RdBu')
        else:
            pl.plot(xcen,ycen,'bo',label='good data')
            pl.plot(xcen[mask],ycen[mask],'ro',label='bad data')
            pl.legend()

    return (xcen,ycen), mask

# def tpf2pix(fname,index=3,rec_array_index=1,verbose=True):
#     '''
#     fname: str, filename
#     rec_array_index: int,
#     index: int, [3,4,5]
#     '''
#     rec_array = fits.open(fname)
#     hdr  = rec_array[rec_array_index].header
#     data = rec_array[rec_array_index].data
#
#     if verbose:
#         objname = hdr['OBJECT']
#         print('Analyzing {}...\n'.format(objname))
#     bjdref = hdr['BJDREFI']
#     start = bjdref+hdr['TSTART']
#     stop  = bjdref+hdr['TSTOP']
#     ndata = len(data)
#
#     times = np.linspace(start,stop,ndata)
#
#     shape = data[0][index].shape
#     h,w = shape[0], shape[1]
#
#     fluxes = np.zeros((len(data),h,w))
#
#     for i in range(len(data)):
#         fluxes[i,:,:] = data[i][index]
#     return times, fluxes
#
#
# def pix2lc(times,fluxes,aper_rad,aper_shape='round',
#            cutoff_limit=1.0,centroid_shift=1,method='abs_distance'):
#     #get centroids
#     centroids, centroid_mask = get_centroids(fluxes, centroid_shift,
#                                         check_outliers=True,method=method)
#     xcen, ycen = np.array(centroids[0]), np.array(centroids[1])
#     #import pdb; pdb.set_trace()
#     x,y = xcen[~centroid_mask],ycen[~centroid_mask]
#
#     aperture_mask = make_mask(fluxes,cutoff_limit,shape=aper_shape)
#
#     flux = []
#     time = np.copy(times)
#
#
#     for n,img in enumerate(fluxes):
#         img = img[aperture_mask]
#         flux.append(np.nansum(img,axis=0))
#     #remove outliers
#     flux = np.array(flux)
#     t,f = time[~centroid_mask], flux[~centroid_mask]
#     df=pd.DataFrame(np.c_[t,f,x,y],index=range(len(t)))
#     df.columns = ['time','flux','Xcen','Ycen']
#     return df


def plot_aper_mask(fluxes,rad,aper_shape,contrast=0.1,epic=None):
    stacked = np.nansum(fluxes,axis=0)
    mask = make_mask(fluxes, rad = rad, shape=aper_shape)
    x,y = centroid_com(stacked)

    fig, ax = pl.subplots(1,1)
    interval = ZScaleInterval(contrast=contrast)
    zmin,zmax = interval.get_limits(stacked)

    cb = ax.imshow(stacked, origin='bottom', interpolation=None)
    ax.imshow(mask, alpha=0.3)
    if aper_shape == 'round':
        circ = pl.Circle((x,y),rad,color='w',alpha=0.2,lw=5,label='r={}'.format(rad))
        ax.add_artist(circ)

    ax.plot(x,y,'r+',ms=20, lw=10,label='centroid')
    pl.colorbar(cb, ax=ax)
    pl.xlabel('X')
    pl.ylabel('Y')
    pl.legend()
    if epic is not None:
        pl.title(epic)
    pl.show()

    return fig

def make_mask(fluxes, rad=None, cutoff_limit=1.0, shape='round', epic=None, showfig=False):
    '''
    create mask given 3D (t x h x w) fluxes
    '''
    #stack images
    stacked = np.nansum(fluxes,axis=0)

    #get centroid of stacked image
    x,y = centroid_com(stacked)
    xx,yy = int(round(x)), int(round(y))

    mask = np.zeros_like(stacked)

    if shape == 'square':
        if rad is not None:
            if stacked.size > (2*rad)**2:
                assert ValueError('mask size bigger than tpf cutout')
            for i in np.arange(-rad,rad,1):
                for j in np.arange(-rad,rad,1):
                    mask[xx+i,yy+j] = 1
        else:
            raise ValueError('rad: aperture radius not defined')

    elif shape == 'irregular':
        mask = find_aperture(fluxes, cutoff_limit=1.0)

    elif shape == 'round' or shape == 'cirle':
        #default: round
        if rad is not None:
            if stacked.size > (2*rad)**2:
                assert ValueError('mask size bigger than tpf cutout')
            for i in range(stacked.shape[0]):
                for j in range(stacked.shape[1]):
                    if np.sqrt((xx-i)**2+(yy-j)**2)<=rad:
                        mask[i,j] = 1
        else:
            raise ValueError('rad: aperture radius not defined')
    else:
        # no mask; use all pixels
        mask = np.ones_like(stacked)

    if showfig:
        fig = plot_aper_mask(fluxes,rad,aper_shape=shape,contrast=0.1,epic=epic)

    return np.array(mask,dtype=bool)

def tpf2lc(fname, radii, aper_shape='round', outlier_sigma=5,
        flat_window=301, corr_window=51, cutoff_limit=1.0, polyorder=4,
        break_tolerance=5,save_as_tpf=False, verbose=False, outdir='reduced'):
    '''
    do aperture photometry with multiple apertures and a mask
    '''
    print('\nAperture photometry with r={} and {} mask...\n'.format(radii,aper_shape))
    if verbose:
        print("""sigma cut for outliers: {}\nwindow length (flatten): {}\nwindow length (sff): {}\ncutoff limit (if mask=irregular): {}\n
            """.format(outlier_sigma,flat_window,corr_window,cutoff_limit))
    hdr = fits.getheader(fname)
    hdulist = fits.open(fname)
    tpf = KeplerTargetPixelFile(fname, quality_bitmask='hardest')
    assert str(tpf.keplerid) in hdulist.filename()

    flux_per_r = {}
    for r in radii:
        mask = make_mask(tpf.flux, rad=r, shape=aper_shape, epic=tpf.keplerid)
        lc = tpf.to_lightcurve(aperture_mask=mask);
        lc2 = lc.remove_nans().remove_outliers(sigma=outlier_sigma)
        flat_lc2, trend = lc2.flatten(window_length=flat_window,
                                    polyorder=polyorder,
                                    break_tolerance=break_tolerance,
                                    return_trend=True)
        corr_lc = flat_lc2.correct(method='sff',windows=corr_window)

        flux_per_r[r]=(corr_lc.time,corr_lc.flux,corr_lc.flux_err)

    if save_as_tpf:
        hdr['cdpp'] = corr_lc.cdpp()
        #append to hdulist photometry of each aperture and save
        for num,r in enumerate(flux_per_r):
            comment_num = 'COMMENT{}'.format(num)
            aper_name = '{}_APER{}'.format(aper_shape,num)
            hdr[comment_num] = 'rad: {}, shape: {}'.format(r,aper_shape)

            tab = table.Table(flux_per_r[r], names=['time','flux','flux_err'])
            bintab=fits.BinTableHDU(tab,name=aper_name,header=hdr)
            #append to original
            hdulist.append(bintab)

        #make hdu for mask
        hdu=fits.hdu.ImageHDU(np.array(mask,dtype=float), name='APERTURE', header=hdr) #problem with bool
        #replace aperture
        hdulist[2] = hdu

        #save fits
        fname_new = os.path.join(outdir,fname.split('/')[-1].split('-')[0][4:]+'_'+aper_shape+'.fits')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        hdulist.writeto(fname_new)
        print('Saved: {}\n'.format(fname_new))

    return flux_per_r, mask


def read_tpf(fname,index,return_hdr=True):
    '''
    fname: str, filename
    index: int, hdulist index
            [0,1,2] = primary, target table, aperture mask
            [4,...] = photometry using specified aperture
    return_hdr: bool
    '''
    hdulist = fits.open(fname)
    if index == 0: #primary
        data = hdulist[index].data
        hdr = hdulist[index].header
        if return_hdr:
            return data, hdr
        else:
            data
    elif index == 1: #target tables
        data = hdulist[index].data
        hdr = hdulist[index].header
        if return_hdr:
            return data, hdr
        else:
            data
    elif index == 2: #aperture mask
        data = hdulist[index].data
        hdr = hdulist[index].header
        if return_hdr:
            return data, hdr
        else:
            data
    else:
        df=table.Table(hdulist[index].data).to_pandas()
        hdr = hdulist[index].header
        if return_hdr:
            return df, hdr
        else:
            df

#---------------------------STATS---------------------------#


def noise_statistic(t, f, timescale=0.25, verbose=False):
    '''
    c.f. lightkurve.cdpp()
    '''
    nchunks = int((t[-1]-t[0])/timescale)+1
    idx = [(t > t[0] + n * timescale) & (t < t[0] + (n + 1) * timescale) for n in range(nchunks)]
    chunks = [f[ix] for ix in idx if ix.sum() > 1]

    cdpp = np.std([np.nanmedian(ch) for ch in chunks])

    if verbose:
        print('cdpp = {:.4f}'.format(cdpp))

    return cdpp
