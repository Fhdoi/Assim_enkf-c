
#===============================================================================
#
# Functions download and interpolate amsr2 6.25km product
# from University of Bremen
# The uncertainty is calculated following Spreen et al. (2008), JGR
#
# Author: Keguang Wang, keguang.wang@met.no
# 25/02/2020
#
#===============================================================================

import os
import sys
import numpy as np
import subprocess as sp
import numpy as np

def amsr2_download(date, wdir):
    web0 = 'https://seaice.uni-bremen.de/data/amsr2/asi_daygrid_swath/n6250/'
    web  = web0 + 'netcdf/' + date[:4] + '/'
    fil  = 'asi-AMSR2-n6250-{}-v5.4.nc'.format(date)
    f1   = os.path.join(wdir, 'amsr2_{}.nc'.format(date))

    if not os.path.exists(f1):
       try:
          print('download AMSR2 ice concentration ...')
          output = sp.check_call(['wget', web + fil])
          if output == 0:
             sp.call(['mv', fil, f1])
             return True
          else:
             return False
       except sp.CalledProcessError:
           print('AMSR2 SIC data not available on ' + date + ' !!!')
           return False
    
    print("AMSR2 data already exists for {}".format(date))
    return True  # return True if it file was already downloaded
       

def amsr2_nearestneighbour(X0,Y0,var0,X,Y,mask,geo=False,nodata=-999):

    print('  interpolating data with nearest-neighbour method ...')

    # get new coordiante dimensions
    ny, nx = X.shape

    #define the origin in old coordinate, and steps
    x0, y0 = X0[0,0], Y0[0,0]
    dx = X0[0,1] - X0[0,0]
    dy = Y0[1,0] - Y0[0,0]
    dd = dx * dy
    #print(dx, dy)

    dimv = len(var0.shape)
    if dimv == 2:
       ny0, nx0 = var0.shape
       var = np.zeros((ny,nx)).astype(np.float32) + np.nan
    elif dimv == 3:
       nt, ny0, nx0 = var0.shape
       var = np.zeros((nt,ny,nx)).astype(np.float32) + np.nan
    elif dimv == 4:
       nt, nz, ny0, nx0 = var0.shape
       var = np.zeros((nt,nz,ny,nx)).astype(np.float32) + np.nan

    # geographical
    #if geo:
    #   X[X < 0] = X[X < 0] + 360.

    # interpolation
    i1 = np.array(np.floor((X-x0) / dx + 0.5)).astype(int)
    i2 = np.array(np.floor((Y-y0) / dy + 0.5)).astype(int)
    #i1[i1 > nx] = nx
    #i2[i2 > ny] = ny

    if dimv == 2:
       var = var0[i2,i1]
       #var[mask == 0] = nodata
    elif dimv == 3:
       var = var0[:,i2,i1]
       for k in range(nt):
           aa = var[k,:,:]
           #aa[mask == 0] = nodata
           var[k,:,:] = aa
    elif dimv == 4:
       var = var0[:,:,i2,i1]
       for k in range(nt):
           for i3 in range(nz):
               aa = var[k,i3,:,:]
               #aa[mask == 0] = nodata
               var[k,i3,:,:] = aa

    return var

def amsr2_err(sic0,mask,nodata=-999):

    print('  calculating uncertainties for amsr2 data ...')

    sic = sic0 + 0.0
    sic[sic == nodata] = 0
    sic = sic * 0.01

    # estimate uncertainty (see Spreen et al., 2008)
    Psw,   Epsw    = 82.0, 4.0
    Psi,   Epsi	   = 10.0, 4.0
    Tauw,  Etauw   = 0.27, 0.1
    Taui,  Etaui   = 0.14, 0.035
    #d3, d2, d1, d0 = 1.64e-5, -1.6e-3, 1.92e-2, 0.971
    d3, d2, d1, d0 = 5.587e-06, -5.218e-04, -1.226e-02, 1.116

    Ps   = sic * Psi   + (1-sic) * Psw
    Tau  = sic * Taui  + (1-sic) * Tauw
    Etau = sic * Etaui + (1-sic) * Etauw
    ac   = 1.1*np.exp(-2*Tau) - 0.11*np.exp(-Tau)
    P    = Ps * ac

    Ep2 = (Ps*Etau*(0.11*np.exp(-Tau)-2.2*np.exp(-2*Tau)))**2 + \
          (ac*(1-sic)*Epsw)**2 + (ac*sic*Epsi)**2
    err = np.abs(3*d3*P**2 + 2*d2*P + d1) * np.sqrt(Ep2) * 100
    #err[mask == 0] = nodata

    return err

def amsr2_writedim(f,Lon,Lat,dt,ensemble=[]):

    print('  writing dimensions ...')

    f.createDimension('time',None)
    time = f.createVariable('time','f8',('time',))
    time.long_name = 'time'
    time.units = 'days since 1970-01-01 00:00:00'
    time[:] = dt #dt.astype(np.float64)

    ne = len(ensemble)
    if ne > 0:
       f.createDimension('number',ne)
       ensid = f.createVariable('number','i4',('number',))
       ensid[:] = ensemble

    if len(Lon.shape) == 2:
       ny, nx = Lon.shape
       f.createDimension('lon',nx)
       f.createDimension('lat',ny)
       lonid = f.createVariable('longitude','f',('lat','lon'))
       latid = f.createVariable('latitude','f',('lat','lon'))
    elif len(Lon.shape) == 1:
       ny, = Lat.shape
       nx, = Lon.shape
       f.createDimension('lon',nx)
       f.createDimension('lat',ny)
       lonid = f.createVariable('longitude','f',('lon',))
       latid = f.createVariable('latitude','f',('lat',))
    else:
       print('geographical dimenion definition no supported !!!')
       sys.exit()

    lonid.units = 'degree East'
    lonid[:] = Lon.astype(np.float32)

    latid.units = 'degree East'
    latid[:] = Lat.astype(np.float32)

def amsr2_writevar(f,var,var_name,unit,long_name,ks=-1,nodata=-999):

    print('  saving ' + var_name + str(var.shape) + ' ...')

    try:
       varid = f.variables[var_name]
       print('    shape of old data:', varid.shape)
       print('    shape of new data:', var.shape)

       kn = var.shape[0]
       varid[ks:ks+kn,...] = var.astype(np.float32)
    except:
       if len(var.shape) <= 3:
          varid = f.createVariable(var_name,'f',('time','lat','lon'))
       elif len(var.shape) == 4:
          varid = f.createVariable(var_name,'f',('time','number','lat','lon'))

       if len(var.shape) == 2:
          ny, nx = var.shape
          var1 = np.zeros([1,ny,nx])
          var1[0,:,:] = var
          varid[:] = var1.astype(np.float32)
       else:
          varid[:] = var.astype(np.float32)
       varid.long_name = long_name
       varid.units = unit
       varid.time = 'time'
       varid.missing_value = nodata

    f.sync()
