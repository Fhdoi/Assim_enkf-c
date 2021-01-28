from pyproj import CRS, Proj, transform
import pyproj
import netCDF4 as nc
import numpy.matlib as npm
import numpy as np
import datetime
import os
import sys
sys.path.insert(0,'/home/sindremf/PHD2/Work/Assim_enkf-c')
import python_tools.pyromstools_cice_amsr2 as damsr
from pyresample.geometry import GridDefinition
from pyresample import geometry, image, SwathDefinition
import xarray as xr
from python_tools.enkf_c_toolbox_parallel import cmd


def download_amsr_Sindre(date, wdir, enkf_c_dir, add_barents = True):
    # date = datetime.datetime(2018,1,1)
    # wdir = '/home/sindremf/PHD2/Work/Test_assimiation/'
    # enkf_c_dir = '/home/sindremf/PHD2/Work/Assim_enkf-c/'
    

    grid_file = enkf_c_dir+'conf/new_grid_ice.nc'

    # Download data with Keguangfunction    
    damsr.amsr2_download(date.strftime('%Y%m%d'), wdir)

    handle = nc.Dataset(wdir+'amsr2_'+date.strftime('%Y%m%d')+'.nc','r')


    X0 = handle['x'][:]
    Y0 = handle['y'][:]
    X02 = npm.repmat(X0, len(Y0), 1)
    Y02 = npm.repmat(Y0,len(X0),1).transpose()
    sic = handle['z'][:]

    # Get lon/lat from projection
    P = Proj('epsg:3411')
    lon2, lat2 = P(X02,Y02,inverse=True)

    
    # Calculate
    sic2 = sic * 0.01

    # estimate uncertainty (see Spreen et al., 2008)
    Psw,   Epsw    = 82.0, 4.0
    Psi,   Epsi	   = 10.0, 4.0
    Tauw,  Etauw   = 0.27, 0.1
    Taui,  Etaui   = 0.14, 0.035
    #d3, d2, d1, d0 = 1.64e-5, -1.6e-3, 1.92e-2, 0.971
    d3, d2, d1, d0 = 5.587e-06, -5.218e-04, -1.226e-02, 1.116

    Ps   = sic2 * Psi   + (1-sic2) * Psw
    Tau  = sic2 * Taui  + (1-sic2) * Tauw
    Etau = sic2 * Etaui + (1-sic2) * Etauw
    ac   = 1.1*np.exp(-2*Tau) - 0.11*np.exp(-Tau)
    P    = Ps * ac

    Ep2 = (Ps*Etau*(0.11*np.exp(-Tau)-2.2*np.exp(-2*Tau)))**2 + \
          (ac*(1-sic2)*Epsw)**2 + (ac*sic2*Epsi)**2
    err2 = np.abs(3*d3*P**2 + 2*d2*P + d1) * np.sqrt(Ep2) * 100


    file = wdir+'AMSR2Full-'+date.strftime('%Y%m%d')+'.nc'

    # Scale to decimal
    siconc = sic.data[:]/100
    err3 = err2.data[:]/100

    # Set nan values to -1
    err3[np.isnan(siconc)] = -1
    siconc[np.isnan(siconc)] = -1

    # Create a new netcdf file

    try:
        os.remove(file)
    except:
        pass
        
    ds = nc.Dataset(file, 'w', format='NETCDF4')

    time = ds.createDimension('time', None)
    times = ds.createVariable('time', 'f4', ('time',))
    times[:] = nc.date2num(date, units='days since 1990-01-01',
                             calendar='gregorian')
    times.units = 'days since 1990-01-01'
    times.calendar='gregorian'


    dx = ds.createDimension('dx', len(X0))
    dy = ds.createDimension('dy', len(Y0))



    dxs = ds.createVariable('dx', 'f4', ('dx',))
    dys = ds.createVariable('dy', 'f4', ('dy',))

    #print(Nens)
    #print(np.arange(0, Nens, 1.0))
    dxs[:] = X0
    dys[:] = Y0

    lon = ds.createVariable('lon', 'f4', ('dy', 'dx',))
    lat = ds.createVariable('lat', 'f4', ('dy', 'dx',))

    lon[:] = lon2
    lat[:] = lat2
    err = ds.createVariable('error_std', 'f4', ('time', 'dy', 'dx',))
    iceconc = ds.createVariable('iceconc', 'f4', ('time', 'dy', 'dx',))
    iceconc[0,:,:] = siconc
    err[0,:,:] = err3

    if add_barents:
        
        # Load the barents grid
        handle2 = xr.open_dataset(grid_file)
        lon_mod = handle2['lon']
        lat_mod = handle2['lat']
        mod_grid_def = geometry.GridDefinition(lons=lon_mod, lats=lat_mod)
        
        dx2 = ds.createDimension('dx2', lon_mod.shape[0])
        dy2 = ds.createDimension('dy2', lon_mod.shape[1])

        #dxs2 = ds.createVariable('dx2', 'f4', ('dx',))
        #dys2 = ds.createVariable('dy2', 'f4', ('dy',))
        
        lon4 = ds.createVariable('lon_bar', 'f4', ('dx2', 'dy2',))
        lat4 = ds.createVariable('lat_bar', 'f4', ('dx2', 'dy2',))

        lon4[:] = lon_mod
        lat4[:] = lat_mod

        obs_grid_def = geometry.GridDefinition(lons=lon2, lats=lat2)
        
        # Convert sicconc to barents
        obs_container = image.ImageContainerNearest(siconc, 
                    obs_grid_def, radius_of_influence=200000000)
        obs_modelgrid = obs_container.resample(mod_grid_def)
        siconc_bar = obs_modelgrid.image_data
        
        # Convert error to barents
        obs_container = image.ImageContainerNearest(err3, 
                    obs_grid_def, radius_of_influence=200000000)
        obs_modelgrid = obs_container.resample(mod_grid_def)
        err3_bar = obs_modelgrid.image_data
        
        bar = ds.createVariable('siconc_bar', 'f4', ('dx2', 'dy2',))
        Err_bar = ds.createVariable('error_bar', 'f4', ('dx2', 'dy2',))
        
        bar[:] = siconc_bar
        Err_bar[:] = err3_bar
        
        ds.close()
        handle2.close()

        handle.close()

    os.remove(wdir+'amsr2_'+date.strftime('%Y%m%d')+'.nc')
    os.rename(file,wdir+'amsr2_'+date.strftime('%Y%m%d')+'.nc')

def prep_osisaf_obs(date,obs_dir,Assim_dir):
    #osisaf_pre = 'ice_conc_nh_ease-125_multi_'
    osisaf_pre = 'ice_conc_nh_polstere-100_multi_'
    osisaf_post = '1200.nc'
    smnd = str(date.month) if date.month > 9 else '0'+str(date.month)
    sday = str(date.day) if date.day > 9 else '0'+str(date.day)

    obs_file = obs_dir+str(date.year)+'/'+smnd+'/'+osisaf_pre+str(date.year)+smnd+sday+osisaf_post

    file_out = Assim_dir +'/obs/OSISAF/this_day.nc'
    
    if os.path.exists(obs_file):
        # Export concentraion and uncertainty
        cmd('ncks -O -v ice_conc,total_uncertainty '+obs_file+' temp_osisaf1.nc')
        # Rename uncertainty to that read by enkf-c
        cmd('ncrename -v total_uncertainty,error_std temp_osisaf1.nc')
        # Change dimension from percent to decimal concentration
        cmd('ncap2 -O -s "ice_conc=ice_conc/100" temp_osisaf1.nc temp_osisaf2.nc')
        cmd('ncap2 -O -s "error_std=error_std/100" temp_osisaf2.nc '+file_out)

        # Fix nan values 
        cmd('ncatted -O -h -a _FillValue,error_std,m,f,-1 '+file_out)
        cmd('ncatted -O -h -a _FillValue,ice_conc,m,f,-1 '+file_out)

        # delete temporary files
        cmd('rm temp_osisaf1.nc temp_osisaf2.nc')
    else:
        print('Osisaf file does not exist, please download the file')
        try:
            os.remove(file_out)
        except:
            pass
