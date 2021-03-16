#Script for downloading observations from assimlation
import datetime
import xarray as xr
from pyresample.geometry import GridDefinition
from pyresample import geometry, image, SwathDefinition
import numpy as np
from numpy.matlib import repmat
import matplotlib.pyplot as plt
from pyproj import Proj
import subprocess


def download_observations(date_list=None,obstype_list=None,obs_dir='./',custom_file=None,keep_both=False):
# Usage
# date_list, should be a list of datetime objects f.eks: date_list = [base + datetime.timedelta(days=x) for x in range(0,365,7)]
# obstype_list should be a list of the observations that should be downloaded, f.eks obs = ['AMSR','SSMIS','CRYO','SMOS','MUR']
# obs_dir should be the directory where observations are stored: f.eks: obs_dir = '/home/sindre/PHD2/Work/Observations/'
# This will create directories obs_dir + 'SMOS' and so on.
# custom_file should be used when you want the observations to be interpolated to a custom grid. For MUR this is obligatory as the 
# original size is so big that it is too much to store generally. THe grid file should contain a 2d variable called lon and one called
# lat which holdes the grid coordinates you want to transform you model into.
# If a custom file is defined keep_both=True means that both the original file and the custom grid file is kept, if false only the 
# custom file is kept.
    for date in date_list:
        if 'AMSR' in obstype_list:
            # Download AMSR
            print('Downloading AMSR for '+date.strftime('%Y%m%d'))
            download_AMSR(date=date,obs_dir=obs_dir,custom_grid_file=custom_file, keep_both=keep_both)
        if 'CRYO' in obstype_list:
            # Download CRYO
            print('Downloading CRYO for '+date.strftime('%Y%m%d'))
            download_CRYO(date=date,obs_dir=obs_dir,custom_grid_file=custom_file, keep_both=keep_both)
        if 'SMOS' in obstype_list:
            # Download SMOS
            print('Downloading SMOS for '+date.strftime('%Y%m%d'))
            download_SMOS(date=date,obs_dir=obs_dir,custom_grid_file=custom_file, keep_both=keep_both)
        if 'MUR' in obstype_list:
            # Download MUR
            print('Downloading MUR for '+date.strftime('%Y%m%d'))
            download_MUR(date=date,obs_dir=obs_dir,custom_grid_file=custom_file)
        if 'SSMIS' in obstype_list:
            # Download MUR
            print('Downloading SSMIS for '+date.strftime('%Y%m%d'))
            download_SSMIS(date=date,obs_dir=obs_dir,custom_grid_file=custom_file, keep_both=keep_both)
            
def cmd(command):
    """Function runs provided command in the system shell.

    Args:
        command (string) : Command to be executed in shell
    Returns:
        result (integer) : Shell returned status of command
    """
    print("> " + command)
    result = subprocess.call(command, shell = True)

    if result != 0:
        print("Command failed: %d" % result)
    else:
        return result

def create_custom_grid(ds,filename,var_name,dl,custom_grid):
    cgrid = xr.open_dataset(custom_grid)
    clon = cgrid['lon']
    clat = cgrid['lat']
    cgrid_def = geometry.GridDefinition(lons=clon, lats=clat)
    
    ogrid_def = geometry.GridDefinition(lons=ds['lon'], lats=ds['lat'])
    
    # Convert data to custom grid
    obs_container = image.ImageContainerNearest(ds[var_name].data[0,:,:], 
                        ogrid_def, radius_of_influence=20000000)
    obs_modelgrid = obs_container.resample(cgrid_def)
    data_out = obs_modelgrid.image_data
    
    obs_container = image.ImageContainerNearest(ds['error_std'].data[0,:,:], 
                        ogrid_def, radius_of_influence=20000000)
    obs_modelgrid = obs_container.resample(cgrid_def)
    err_out = obs_modelgrid.image_data
    
    ds = xr.Dataset(
            {var_name: (('time','x', 'y'), np.expand_dims(data_out, axis=0)),
            "error_std": (('time','x', 'y'), np.expand_dims(err_out, axis=0)),
            "lon": (('x', 'y'), clon),
            "lat": (('x', 'y'), clat)},
            coords={
                'time': dl[0:1],
            },
            )
    ds.to_netcdf(filename)
    





def download_AMSR(date,obs_dir,custom_grid_file=None, keep_both=True):
    # Could also simply be made to use a custom grid, will consider this, only needs the
    # denne må nok også tas inn som input
    dl = [date]
    amsr_dir = obs_dir+'/AMSR/'
    amsr_pre = 'https://seaice.uni-bremen.de/data/amsr2/asi_daygrid_swath/n6250/netcdf/'
    amsr_mid = '/asi-AMSR2-n6250-'

    #Download the file
    res = cmd('wget -P '+amsr_dir+' '+amsr_pre+date.strftime('%Y')+amsr_mid+date.strftime('%Y%m%d')+'-v5.4.nc')

    # If file exists
    if res == 0:
        amsr_file = amsr_dir + amsr_mid+date.strftime('%Y%m%d')+'-v5.4.nc'


        # Get the coordinates
        handle = xr.open_dataset(amsr_file)

        X0 = handle['x'][:]
        Y0 = handle['y'][:]
        X02 = repmat(X0, len(Y0), 1)
        Y02 = repmat(Y0,len(X0),1).transpose()
        sic = handle['z'][:]

        # Get lon/lat from projection
        P = Proj('epsg:3411')
        lon2, lat2 = P(X02,Y02,inverse=True)

        # Adjust sic
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

        # Scale to decimal
        siconc = sic.data[:]/100
        err3 = err2.data[:]/100

        # Set nan values to -1
        err3[np.isnan(siconc)] = 10
        siconc[np.isnan(siconc)] = -1

        # Change water uncertainty to 0.05, it is way to big
        err3[siconc == 0] = 0.05   
        
        # save the file
        amsr_file2 = amsr_dir+'AMSR2_'+date.strftime('%Y%m%d')+'.nc'
        ds = xr.Dataset(
                {"ice_conc": (('time','x', 'y'), np.expand_dims(siconc, axis=0)),
                "error_std": (('time','x', 'y'), np.expand_dims(err3, axis=0)),
                "lon": (('x', 'y'), lon2),
                "lat": (('x', 'y'), lat2)},
                coords={
                    'time': dl[0:1],
                },
                )
        ds.to_netcdf(amsr_file2)
        
        # Delete original download
        handle.close()
        cmd('rm '+amsr_file)

        if custom_grid_file:
            amsr_file3 = amsr_dir+'AMSR2_'+date.strftime('%Y%m%d')+'b.nc'
            create_custom_grid(ds=ds,filename=amsr_file3,var_name='ice_conc',dl=dl,custom_grid=custom_grid_file)
            if keep_both:
                pass
            else:
                cmd('rm '+amsr_file2)




def download_CRYO(date,obs_dir,custom_grid_file=None, keep_both=True):
    # To download the cryosat data a user to earthdata login is requires
    # see link for more info: https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+cURL+And+Wget
    dl = [date]
    cryo_usr = '--load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies'
    cryo_dir = obs_dir+'/CRYO/'
    cryo_txt = 'https://n5eil01u.ecs.nsidc.org/ICEBRIDGE/RDEFT4.001/'
    
    date2 = date + datetime.timedelta(days=30)
    cryo_file = cryo_dir + 'RDEFT4_'+date2.strftime('%Y%m%d')+'.nc'
    res = cmd('wget -P '+cryo_dir+' '+cryo_usr+' '+cryo_txt+date.strftime('%Y')+'.'+date.strftime('%m')+'.'+date.strftime('%d')+'/RDEFT4_'+date2.strftime('%Y%m%d')+'.nc')

    # If file exists
    if res == 0:
        
        cryo_out = cryo_dir+'CRYO_'+date.strftime('%Y%m%d')+'.nc'
        cmd('ncks -O -v lon,lat,sea_ice_thickness '+cryo_file+' '+cryo_out)
        cmd('rm '+cryo_file)
        cmd('ncrename -O -v sea_ice_thickness,sit '+cryo_out+' '+cryo_out)
        cmd('ncap2 -O -s "where(lon > 180.0) lon = lon-360.0" '+cryo_out+' '+cryo_out)
        cmd('ncap2 -O -s "where(sit < 0.0) sit = -1" '+cryo_out+' '+cryo_out)
        
        
        
        if custom_grid_file:
            ds = xr.open_dataset(cryo_out)
            cryo_out2 = cryo_dir+'CRYO_'+date.strftime('%Y%m%d')+'b.nc'
            ds2 = xr.Dataset(
                {"sit": (('time','x', 'y'), np.expand_dims(ds['sit'].data[:],axis=0)),
                "error_std": (('time','x', 'y'), np.expand_dims(0.1*ds['sit'].data[:],axis=0)),
                "lon": (('x', 'y'), ds['lon'].data[:]),
                "lat": (('x', 'y'), ds['lat'].data[:])},
                coords={
                    'time': dl[0:1],
                },
                )
            
            
            create_custom_grid(ds=ds2,filename=cryo_out2,var_name='sit',dl=dl,custom_grid=custom_grid_file)
            if keep_both:
                pass
            else:
                cmd('rm '+cryo_out)

def download_SMOS(date,obs_dir,custom_grid_file=None, keep_both=True):
    dl = [date]
    smos_pretext = 'http://icdc.cen.uni-hamburg.de/thredds/fileServer/ftpthredds/smos_sea_ice_thickness/v3.2/nh/'
    smos_mid = '/SMOS_Icethickness_v3.2_north_'
    smos_dir = obs_dir+'/SMOS/'
    smos_file = smos_dir + smos_mid + date.strftime('%Y%m%d') +'.nc'
    smos_text = smos_pretext + date.strftime('%Y') +'/'+date.strftime('%m') + smos_mid + date.strftime('%Y%m%d') +'.nc'
    res = cmd('wget -P '+smos_dir+' '+smos_text)
    
    # If file exists
    if res == 0:
        
        # Extract the variables I want
        smos_out = smos_dir+'SMOS_'+date.strftime('%Y%m%d')+'.nc'
        cmd('ncks -O -v longitude,latitude,ice_thickness_uncertainty,sea_ice_thickness '+smos_file+' '+smos_out)
        # Delete the original file
        cmd('rm '+smos_file)
        # Rewrite variable names
        cmd('ncrename -O -v ice_thickness_uncertainty,error_std '+smos_out+' '+smos_out)
        cmd('ncrename -O -v sea_ice_thickness,sit '+smos_out+' '+smos_out)
        cmd('ncrename -O -v longitude,lon '+smos_out+' '+smos_out)
        cmd('ncrename -O -v latitude,lat '+smos_out+' '+smos_out)
        cmd('ncap2 -O -s "where(error_std == 0) error_std = 0.01" '+smos_out+' '+smos_out)
        
        if custom_grid_file:
            ds = xr.open_dataset(smos_out)
            smos_out2 = smos_dir+'SMOS_'+date.strftime('%Y%m%d')+'b.nc'
            
            tt = ds['sit'].data[:]
            et = ds['error_std'].data[:]
            
            et[np.isnan(tt)] = 10
            tt[np.isnan(tt)] = -1
            ds2 = xr.Dataset(
                {"sit": (('time','x', 'y'), ds['sit'].data[:]),
                "error_std": (('time','x', 'y'), ds['error_std'].data[:]),
                "lon": (('x', 'y'), ds['lon'].data[:]),
                "lat": (('x', 'y'), ds['lat'].data[:])},
                coords={
                    'time': dl[0:1],
                },
                )
            
            
            create_custom_grid(ds=ds2,filename=smos_out2,var_name='sit',dl=dl,custom_grid=custom_grid_file)
            if keep_both:
                pass
            else:
                cmd('rm '+smos_out)


def download_MUR(date,obs_dir,custom_grid_file=None):
    # Because of size there is no option for storing the full original file
    # Same as for downloading the CRYO data there is a need for a username a earth data login
    # https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+cURL+And+Wget
    # Because the original MUR files are so large the only option currently is to get the observations
    # on the custom grid provided.
    cryo_usr = '--load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies'
    mur_dir = obs_dir+'/MUR_new/'
    mur_postxt = '090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc'
    mur_pretxt = 'https://podaac-opendap.jpl.nasa.gov/opendap/allData/ghrsst/data/GDS2/L4/GLOB/JPL/MUR/v4.1/'
    dl = [date]
    mur_file = mur_dir + date.strftime('%Y%m%d')+mur_postxt
    res = cmd('wget -P '+mur_dir+' '+cryo_usr+' '+mur_pretxt+date.strftime('%Y')+'/'+date.strftime('%j')+'/'+date.strftime('%Y%m%d')+mur_postxt)

    # If file exists
    if res == 0:
        cmd('ncks -O -d lat,13000,17998 '+mur_file+' '+mur_file)
        cmd('ncks -O -v analysed_sst,analysis_error '+mur_file+' '+mur_file)
        # Just make a new file as I need to convert the grid and also make new lon/lat files.
        DS = xr.open_dataset(mur_file)
        murlon = DS['lon']
        murlat = DS['lat']

        murlat2 = np.transpose(repmat(murlat.data,len(murlon),1))
        murlon2 = repmat(murlon.data,len(murlat),1)

        mur_grid_def = geometry.GridDefinition(lons=murlon2, lats=murlat2)

        sst = DS['analysed_sst']
        sst_err = DS['analysis_error']

        #Read custom grid
        DS_bar = xr.open_dataset(custom_grid_file)
        barlon = DS_bar['lon']
        barlat = DS_bar['lat']
        bar_grid_def = geometry.GridDefinition(lons=barlon, lats=barlat)

        # Convert to barents
        obs_container = image.ImageContainerNearest(sst.data[0,:,:], 
                            mur_grid_def, radius_of_influence=2000)
        obs_modelgrid = obs_container.resample(bar_grid_def)
        sst_bar = obs_modelgrid.image_data
        obs_container = image.ImageContainerNearest(sst_err.data[0,:,:], 
                            mur_grid_def, radius_of_influence=2000)
        obs_modelgrid = obs_container.resample(bar_grid_def)
        sst_err2 = obs_modelgrid.image_data
        # Sematics

        # Remove nan values 
        sst_bar = sst_bar -273.15
        sst_bar[np.isnan(sst_bar)] = -3 
        sst_err2[np.isnan(sst_err2)] = 10
        # Write to file
        ds = xr.Dataset(
                {"sst": (('time','x', 'y'), np.expand_dims(sst_bar, axis=0)),
                "error_std": (('time','x', 'y'), np.expand_dims(sst_err2, axis=0)),
                "lon": (('x', 'y'), barlon),
                "lat": (('x', 'y'), barlat)},
                coords={
                    'time': dl[0:1],
                },
                )
        ds.to_netcdf(mur_dir+'MUR_'+date.strftime("%Y%m%d")+'.nc')
        cmd('rm '+mur_file)

def download_SSMIS(date,obs_dir,custom_grid_file=None,keep_both=False):
    # Ta fra osisaf,
    # cmd('wget ftp://osisaf.met.no/archive/ice/conc/2018/01/ice_conc_nh_polstere-100_multi_201801011200.nc')
    dl = [date]
    ssmis_dir = obs_dir+'SSMIS/'
    ssmis_pre= ' ftp://osisaf.met.no/archive/ice/conc/'
    ssmis_mid = '/ice_conc_nh_polstere-100_multi_'
    res = cmd('wget -P ' + ssmis_dir + ssmis_pre+date.strftime('%Y')+'/'+date.strftime('%m')+'/'+ssmis_mid+date.strftime('%Y%m%d')+'1200.nc')

    # If file exists
    if res == 0:
        ssmis_file = ssmis_dir + ssmis_mid+date.strftime('%Y%m%d')+'1200.nc'
        ssmis_out = ssmis_dir+'SSMIS_'+date.strftime('%Y%m%d')+'.nc'
        cmd('ncks -O -v lon,lat,ice_conc,total_uncertainty '+ssmis_file+' '+ssmis_out)
        cmd('rm '+ssmis_file)
        cmd('ncrename -O -v total_uncertainty,error_std '+ssmis_out+' '+ssmis_out)
        cmd('ncap2 -O -s "ice_conc=ice_conc/100" '+ssmis_out+' '+ssmis_out)
        cmd('ncap2 -O -s "error_std=error_std/100" '+ssmis_out+' '+ssmis_out)
        if custom_grid_file:
            ds = xr.open_dataset(ssmis_out)
            ssmis_out2 = ssmis_dir+'SSMIS_'+date.strftime('%Y%m%d')+'b.nc'
            
            tt = ds['ice_conc'].data[:]
            et = ds['error_std'].data[:]
            
            et[np.isnan(tt)] = 10
            tt[np.isnan(tt)] = -1
            ds2 = xr.Dataset(
                {"ice_conc": (('time','x', 'y'), ds['ice_conc'].data[:]),
                "error_std": (('time','x', 'y'), ds['error_std'].data[:]),
                "lon": (('x', 'y'), ds['lon'].data[:]),
                "lat": (('x', 'y'), ds['lat'].data[:])},
                coords={
                    'time': dl[0:1],
                },
                )
            
            
            create_custom_grid(ds=ds2,filename=ssmis_out2,var_name='ice_conc',dl=dl,custom_grid=custom_grid_file)
            if keep_both:
                pass
            else:
                cmd('rm '+ssmis_out)
        
