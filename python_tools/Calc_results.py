import xarray as xr
import numpy as np
import datetime
from os.path import exists
import multiprocessing
from joblib import Parallel, delayed

def read_ens(j,var,file_str, pre_str,date,typed,day=-1,lvl=-10):
    if typed == 'ice':
        file = file_str+str(j).zfill(3)+pre_str+date.strftime("%Y-%m-%d")+'-00000.nc'
    elif typed == 'ocean':
        date2 = date - datetime.timedelta(days=7)
        file = file_str+str(j).zfill(3)+pre_str+date2.strftime("%Y%m%d")+'_'+ date.strftime("%Y%m%d") +'.nc'
    ds2 = xr.open_dataset(file)
    
    #If only one date in the file (ice) not for ocean
    #print(day)
    if day == -1:
        A = ds2[var][:].data
    else:
        # if a specific level is defined
        #print('read')
        #print(lvl)
        if lvl > -2:
            #print(j)
            A = ds2[var][day,lvl,...].data
        else:
            A = ds2[var][day,...].data                
        
    ds2.close()
    return np.squeeze(A)
    

def write_ice_hist(pre_str,file_str,histdir,date,ice_omit,Nens,nproc=11):
    filename = file_str+str(1).zfill(3)+pre_str+date.strftime("%Y-%m-%d")+'-00000.nc'
    
    ds = xr.open_dataset(filename)
    

    for var in ds.variables:
        if var not in ice_omit:
            A = np.zeros(ds[var].shape+(Nens,))
            A = Parallel(n_jobs=nproc)(delayed(read_ens)(j,var,file_str, pre_str,date,typed='ice',day=-1) for j in range(1, Nens+1))
            #for j in range(2, Nens+1):
            #    A[...,j-1] = read_ens(j,var,file_str, pre_str)
            
            A = np.array(A)
            mean = np.mean(A,0)
            std = np.std(A,0)

            if len(A.shape) == 3:
                ds_res = xr.Dataset(
                {var+'_mean': (('x', 'y'), mean),
                 var+'_std': (('x', 'y'), std),
                },
                )
            elif len(A.shape) == 4:
                ds_res = xr.Dataset(
                {var+'_mean': (('z','x', 'y'), mean),
                 var+'_std': (('z','x', 'y'), std),
                },
                )       
            
            if exists(histdir+'Ensemble_stat_hist'+date.strftime("%Y%m%d")+'.nc'):
                ds_res.to_netcdf(histdir+'Ensemble_stat_hist'+date.strftime("%Y%m%d")+'.nc', mode = 'a')
            else:
                ds_res.to_netcdf(histdir+'Ensemble_stat_hist'+date.strftime("%Y%m%d")+'.nc', mode = 'w')
    ds.close()
    
    
    
def write_ocean_hist(pre_str,file_str,histdir,date,ocean_var,Nens,nproc=11): 

    date2 = date - datetime.timedelta(days=7)
    filename = file_str+str(1).zfill(3)+pre_str+date2.strftime("%Y%m%d")+'_'+ date.strftime("%Y%m%d") +'.nc'
    ds = xr.open_dataset(filename)

    for var in ds.variables:
        if var in ocean_var:
            print(var)
            for day in range(1,8):
                date3 = date2 + datetime.timedelta(days=day)
                A = Parallel(n_jobs=nproc)(delayed(read_ens)(j,var,file_str, pre_str,date,'ocean',day=day) for j in range(1, Nens+1))
                
                A = np.array(A)
                dateut = date + datetime.timedelta(days=day)
                mean = np.mean(A,0)
                std = np.std(A,0)
                
                if var == 'u' or var == 'ubar' or var == 'sustr':
                    xvar = 'xu'
                    yvar = 'yu'
                elif var == 'v' or var == 'vbar' or var == 'svstr':
                    xvar = 'xv'
                    yvar = 'yv'
                else:
                    xvar = 'x'
                    yvar = 'y'  
                    
                zvar = 'z'
                if var == 'w':
                    zvar = 'zw'

                if len(A.shape) == 3:
                    ds_res = xr.Dataset(
                    {var+'_mean': ((xvar, yvar), mean),
                     var+'_std': ((xvar, yvar), std),
                    },
                    )
                elif len(A.shape) == 4:
                    ds_res = xr.Dataset(
                    {var+'_mean': ((zvar,xvar, yvar), mean),
                     var+'_std': ((zvar,xvar, yvar), std),
                    },
                    )
                if exists(histdir+'Ensemble_stat_hist'+date3.strftime("%Y%m%d")+'.nc'):
                    ds_res.to_netcdf(histdir+'Ensemble_stat_hist'+date3.strftime("%Y%m%d")+'.nc', mode = 'a')
                else:
                    ds_res.to_netcdf(histdir+'Ensemble_stat_hist'+date3.strftime("%Y%m%d")+'.nc', mode = 'w')
                ds_res.close()
    ds.close()
                
def write_ocean_surface(pre_str,file_str,histdir,date,ocean_var,var_4d,Nens,nproc=11): 

    date2 = date - datetime.timedelta(days=7)
    filename = file_str+str(1).zfill(3)+pre_str+date2.strftime("%Y%m%d")+'_'+ date.strftime("%Y%m%d") +'.nc'
    ds = xr.open_dataset(filename)

    for var in ds.variables:
        
        if var in ocean_var:
            print(var)
            lvl = -10
            if var in var_4d:
                lvl = -1
            #print(lvl)
            for day in range(1,8):
                date3 = date2 + datetime.timedelta(days=day)
                #print(day)
                A = Parallel(n_jobs=nproc)(delayed(read_ens)(j,var,file_str, pre_str,date,'ocean',day=day,lvl=lvl) for j in range(1, Nens+1))
                #print(day)
                A = np.array(A)
                dateut = date + datetime.timedelta(days=day)


                if var == 'u' or var == 'ubar':
                    xvar = 'xu'
                    yvar = 'yu'
                elif var == 'v' or var == 'vbar':
                    xvar = 'xv'
                    yvar = 'yv'
                else:
                    xvar = 'x'
                    yvar = 'y'  

                ds_res = xr.Dataset(
                {var: (('Nens',xvar, yvar), A),
                },
                )
                if exists(histdir+'Ensemble_surface'+date3.strftime("%Y%m%d")+'.nc'):
                    ds_res.to_netcdf(histdir+'Ensemble_surface'+date3.strftime("%Y%m%d")+'.nc', mode = 'a')
                else:
                    ds_res.to_netcdf(histdir+'Ensemble_surface'+date3.strftime("%Y%m%d")+'.nc', mode = 'w')
                ds_res.close()
    ds.close()
              
def write_ice_surface(pre_str,file_str,histdir,date,ice_omit,Nens,nproc=11):
    filename = file_str+str(1).zfill(3)+pre_str+date.strftime("%Y-%m-%d")+'-00000.nc'
    
    ds = xr.open_dataset(filename)
    

    for var in ds.variables:
        if var not in ice_omit:
            print(var)
            A = np.zeros(ds[var].shape+(Nens,))
            A = Parallel(n_jobs=nproc)(delayed(read_ens)(j,var,file_str, pre_str,date,typed='ice',day=-1) for j in range(1, Nens+1))
            
            A = np.array(A)
            ds_res = xr.Dataset(
            {var: (('Nens','x', 'y'),A),
            },
            )
             
            if exists(histdir+'Ensemble_surface'+date.strftime("%Y%m%d")+'.nc'):
                ds_res.to_netcdf(histdir+'Ensemble_surface'+date.strftime("%Y%m%d")+'.nc', mode = 'a')
            else:
                ds_res.to_netcdf(histdir+'Ensemble_surface'+date.strftime("%Y%m%d")+'.nc', mode = 'w')
            print(histdir+'Ensemble_surface'+date.strftime("%Y%m%d")+'.nc')
            ds_res.close()
    ds.close()
    
   
