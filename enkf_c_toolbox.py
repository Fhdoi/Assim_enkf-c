# This toolbox contains functions for prepping, tuning and updating when using the enkf-c
# for assimilation.

# There are still som uncertainties regarding the tuning that needs to be solve when the
# format of the observation files are decided, but only one sst and one ice conc file is 
# used as the present plan is, the tuning should be working well.

import glob
import netCDF4 as nc
import datetime
import numpy as np
import os
import sys
import copy

from scipy.ndimage import uniform_filter

# Rewrite the model output to files recognized by the enkf-c
def Prep_ensemble(ens_date, grid_dir, ens_inn_dir, enkf_c_dir, res_type, EnKF_var,Nens):
    # grid_dir: file containing longitude and latitude of the model
    # ens_inn_dir: folder containning the model output data
    # enkf-c_dir: folder containting the netcdf files used by enkf-c
    # ens_inn_file_dummy: file name shared by all the ensemble members
    ens_count = 0
    
    file_cord_handle = nc.Dataset(grid_dir, mode='r')
    
    ice_halo_cells = False
    # Disse leses fra grid fila!    
    lat_rho = file_cord_handle.variables['lat_rho']
    #lon_rho = file_cord_handle.variables['lon_rho']

    # Write the files that are used to disk
    #glob.glob(folder+file_dummy)
    #for ff in glob.glob(folder+file_dummy):
    #    file1.writelines(ff+'\n')
    prescripts = ['iced.','ocean.']
    syear = str(ens_date.year)
    smnd = str(ens_date.month)
    if ens_date.month < 10:
        smnd = '0' + smnd
    sday = str(ens_date.day)
    if ens_date.day < 10:
        sday = '0' + sday

    # Write the files that are used to disk
    #glob.glob(folder+file_dummy)
    file_ens = open(enkf_c_dir+'files_in_ensemble', "w")

    for pre in prescripts:
        iens2 = 0
        for iens in range(1,Nens+1):
            s1ens = str(iens)
            if iens < 100:
                s1ens = '0'+s1ens
            if iens < 10:
                s1ens = '0'+s1ens
            file = ens_inn_dir+pre+syear+smnd+sday+'_'+s1ens+'.nc'
            print(file)
            #This accounts for the possibilty that some ensemble members have not finished in time
            if os.path.exists(file):
                

                iens2 += 1
                sens = str(iens2)
                if iens2 < 100:
                    sens = '0'+sens
                if iens2 < 10:
                    sens = '0'+sens
                
                file_handle = nc.Dataset(file, mode='r')
                if pre == 'ocean.':
                    # Write the file name to the member files folder
                    # No need to wrote it two times as both files are assumed to exist
                    file_ens.writelines(s1ens+'\n')
                    # The ocean restart files can contain several restart times.
                    file_time = file_handle.variables['ocean_time']
                    target_num = nc.date2num(ens_date, units=file_time.units,
                                 calendar=file_time.calendar)
                    date_index = np.where(np.array(file_time[:]) == target_num)[0]
                
                elif pre == 'iced.':
                    date_index = 0
                    test = file_handle.variables['uvel']
                    # Check the if the sea ice restart files is using haloc cells for boundary conditions
                    if test.shape[0]>lat_rho.shape[0]:
                        ice_halo_cells = True


       



                # read the important variables from the file
                for var in file_handle.variables.keys():
                    #print(var)

                    if var in EnKF_var:

      

                        fn = enkf_c_dir+'ensemble_6565/mem'+sens+'_'+var+'.nc'
                        ds = nc.Dataset(fn, 'w', format='NETCDF4')

                        var_inn = file_handle.variables[var]

                        # Litt usikker på om jeg skal ha time, 
                        # teste om det har noe å si på assimilasjonen
                        time = ds.createDimension('time', None)

                        if len(var_inn.shape) == 4:
            
                            dx = ds.createDimension('dx', var_inn.shape[2])
                            dy = ds.createDimension('dy', var_inn.shape[3])
                            dz = ds.createDimension('dz', var_inn.shape[1])

                            times = ds.createVariable('time', 'f4', ('time',))
                            dxs = ds.createVariable('dx', 'f4', ('dx',))
                            dys = ds.createVariable('dy', 'f4', ('dy',))
                            dzs = ds.createVariable('dz', 'f4', ('dz',))

                            temps = ds.createVariable(var, 'f4', ('time', 'dz','dx', 'dy',))

                            dxs[:] = np.arange(0, var_inn.shape[2], 1.0)
                            dys[:] = np.arange(0, var_inn.shape[3], 1.0)
                            dzs[:] = np.arange(0, var_inn.shape[1], 1.0)
                            
                            if pre == 'ocean.':
                                temps[:,:,:,:] = var_inn[date_index,:,:,:]
                            else:
                                sys.exit('Not implemented for ice!')

                        elif len(var_inn.shape) == 3:
                            

                            times = ds.createVariable('time', 'f4', ('time',))
                            if ice_halo_cells:
                                dx = ds.createDimension('dx', var_inn.shape[1]-2)
                                dy = ds.createDimension('dy', var_inn.shape[2]-2)

                                dxs = ds.createVariable('dx', 'f4', ('dx',))
                                dys = ds.createVariable('dy', 'f4', ('dy',))
                                
                                dxs[:] = np.arange(0, var_inn.shape[1]-2, 1.0)
                                dys[:] = np.arange(0, var_inn.shape[2]-2, 1.0)
                            else:
                                dx = ds.createDimension('dx', var_inn.shape[1])
                                dy = ds.createDimension('dy', var_inn.shape[2])

                                dxs = ds.createVariable('dx', 'f4', ('dx',))
                                dys = ds.createVariable('dy', 'f4', ('dy',))
                                
                                dxs[:] = np.arange(0, var_inn.shape[1], 1.0)
                                dys[:] = np.arange(0, var_inn.shape[2], 1.0)

                            if pre == 'ocean.':
                                temps = ds.createVariable(var, 'f4', ('time', 'dx', 'dy',))
                                temps[:,:,:] = var_inn[date_index,:,:]
                            elif pre == 'iced.':
                                dz = ds.createDimension('dz', var_inn.shape[0])
                                dzs = ds.createVariable('dz', 'f4', ('dz',))
                                temps = ds.createVariable(var, 'f4', ('time','dz', 'dx', 'dy',))
                                dzs[:] = np.arange(0, var_inn.shape[0], 1.0)
                
                                
                                if ice_halo_cells:
                                    temps[0,:,:,:] = var_inn[:,1:-1,1:-1]
                                else: 
                                    temps[0,:,:,:] = var_inn[:,:,:]

                        elif len(var_inn.shape) == 2:
                            time = ds.createDimension('time', None)
                            

                            times = ds.createVariable('time', 'f4', ('time',))
        
                            if ice_halo_cells:
                                dx = ds.createDimension('dx', var_inn.shape[0]-2)
                                dy = ds.createDimension('dy', var_inn.shape[1]-2)

                                dxs = ds.createVariable('dx', 'f4', ('dx',))
                                dys = ds.createVariable('dy', 'f4', ('dy',))

                                dxs[:] = np.arange(0, var_inn.shape[0]-2, 1.0)
                                dys[:] = np.arange(0, var_inn.shape[1]-2, 1.0)
                            else: 
                                dx = ds.createDimension('dx', var_inn.shape[0])
                                dy = ds.createDimension('dy', var_inn.shape[1])

                                dxs = ds.createVariable('dx', 'f4', ('dx',))
                                dys = ds.createVariable('dy', 'f4', ('dy',))

                                dxs[:] = np.arange(0, var_inn.shape[0], 1.0)
                                dys[:] = np.arange(0, var_inn.shape[1], 1.0)

                            temps = ds.createVariable(var, 'f4', ('time', 'dx', 'dy',))

      



                            if pre == 'ocean.':
                                sys.exit('Not implemented for ocean!')
                            elif pre == 'iced.':
                                if ice_halo_cells:
                                    temps[0,:,:] = var_inn[1:-1,1:-1]
                                else: 
                                    
                                    temps[0,:,:] = var_inn[:,:]
                                

                        #lat = ds.createVariable('lat', 'f4', ('dx', 'dy',))
                        #lon = ds.createVariable('lon', 'f4', ('dx', 'dy',))

                        times[:] = nc.date2num(ens_date, units='days since 1990-01-01',
                             calendar='gregorian')
                        times.units = 'days since 1990-01-01'
                        times.calendar='gregorian'

                        if var == 'temp':
                            temps.units = 'degrees C'

                        #lat[:,:] = lat_rho[:,:]
                        #lon[:,:] = lon_rho[:,:]

                        #value[0,:,:] = 3

                        ds.close()

                        if var == 'temp' or var == 'salt':
                        # Write also sst variable, must figure out which depth it is and then calculate the variable to
                        # the surface if possible, Johannes should know more about this

                            #print('Also make an sst parameter')
                            
                            if var == 'temp':
                                var2 = 'sst'
                                fn = enkf_c_dir+'ensemble_6565/mem'+sens+'_sst.nc'
                            if var == 'salt':
                                var2 = 'sss'
                                fn = enkf_c_dir+'ensemble_6565/mem'+sens+'_sss.nc'
                            ds = nc.Dataset(fn, 'w', format='NETCDF4')

                            #print(var)
                            time = ds.createDimension('time', None)
                            dx = ds.createDimension('dx', var_inn.shape[2])
                            dy = ds.createDimension('dy', var_inn.shape[3])


                            times = ds.createVariable('time', 'f4', ('time',))
                            dxs = ds.createVariable('dx', 'f4', ('dx',))
                            dys = ds.createVariable('dy', 'f4', ('dy',))
                            sst = ds.createVariable(var2, 'f4', ('time', 'dx', 'dy',))
                        #    lat = ds.createVariable('lat', 'f4', ('dx', 'dy',))
                        #    lon = ds.createVariable('lon', 'f4', ('dx', 'dy',))
                            if var == 'temp':
                                sst.units = 'degrees C'

                            times[:] = file_time[date_index]
                            times.units = file_time.units

                            dxs[:] = np.arange(0, var_inn.shape[2], 1.0)
                            dys[:] = np.arange(0, var_inn.shape[3], 1.0)
                            sst[0,:,:] = var_inn[0,41,:,:]
                            #lat[:,:] = lat_rho[:,:]
                            #lon[:,:] = lon_rho[:,:]

                            #value[0,:,:] = 3

                            ds.close()

                        if var == 'aicen' or var == 'vicen':
                        # Write also integrated variables for assimilation
                            
                            if var == 'aicen':
                                var2 = 'aice'
                                fn = enkf_c_dir+'ensemble_6565/mem'+sens+'_aice.nc'
                            elif var == 'vicen':
                                var2 = 'vice'
                                fn = enkf_c_dir+'ensemble_6565/mem'+sens+'_vice.nc'
                            ds = nc.Dataset(fn, 'w', format='NETCDF4')

                            time = ds.createDimension('time', None)
                            times = ds.createVariable('time', 'f4', ('time',))



                            if ice_halo_cells:
                                dx = ds.createDimension('dx', var_inn.shape[1]-2)
                                dy = ds.createDimension('dy', var_inn.shape[2]-2)

                                dxs = ds.createVariable('dx', 'f4', ('dx',))
                                dys = ds.createVariable('dy', 'f4', ('dy',))

                                dxs[:] = np.arange(0, var_inn.shape[1]-2, 1.0)
                                dys[:] = np.arange(0, var_inn.shape[2]-2, 1.0)
                            else: 
                                dx = ds.createDimension('dx', var_inn.shape[1])
                                dy = ds.createDimension('dy', var_inn.shape[2])

                                dxs = ds.createVariable('dx', 'f4', ('dx',))
                                dys = ds.createVariable('dy', 'f4', ('dy',))

                                dxs[:] = np.arange(0, var_inn.shape[1], 1.0)
                                dys[:] = np.arange(0, var_inn.shape[2], 1.0)

                            sst = ds.createVariable(var2, 'f4', ('time', 'dx', 'dy',))
                            
                            if ice_halo_cells:
                                sst[0,:,:] = np.sum(var_inn[:,1:-1,1:-1],axis=0)
                            else:
                                sst[0,:,:] = np.sum(var_inn[:,:,:],axis=0)
                            #lat[:,:] = lat_rho[:,:]
                            #lon[:,:] = lon_rho[:,:]

                            #value[0,:,:] = 3

                            times[:] = nc.date2num(ens_date, units='days since 1990-01-01',
                             calendar='gregorian')
                            times.units = 'days since 1990-01-01'
                            times.calendar='gregorian'

                            ds.close()

                        if var == 'aicen':
                        # Write smoothed variable for low resolution assimilation
                            
                            fn = enkf_c_dir+'ensemble_6565/mem'+sens+'_aiceosi.nc'
                            ds = nc.Dataset(fn, 'w', format='NETCDF4')

                            time = ds.createDimension('time', None)
                            times = ds.createVariable('time', 'f4', ('time',))



                            if ice_halo_cells:
                                dx = ds.createDimension('dx', var_inn.shape[1]-2)
                                dy = ds.createDimension('dy', var_inn.shape[2]-2)

                                dxs = ds.createVariable('dx', 'f4', ('dx',))
                                dys = ds.createVariable('dy', 'f4', ('dy',))

                                dxs[:] = np.arange(0, var_inn.shape[1]-2, 1.0)
                                dys[:] = np.arange(0, var_inn.shape[2]-2, 1.0)
                            else: 
                                dx = ds.createDimension('dx', var_inn.shape[1])
                                dy = ds.createDimension('dy', var_inn.shape[2])

                                dxs = ds.createVariable('dx', 'f4', ('dx',))
                                dys = ds.createVariable('dy', 'f4', ('dy',))

                                dxs[:] = np.arange(0, var_inn.shape[1], 1.0)
                                dys[:] = np.arange(0, var_inn.shape[2], 1.0)

                            sst = ds.createVariable('aiceosi', 'f4', ('time', 'dx', 'dy',))
                            
                            if ice_halo_cells:
                                sst[0,:,:] = uniform_filter(np.sum(var_inn[:,1:-1,1:-1],axis=0), size=16, mode='constant')
                            else:
                                sst[0,:,:] = uniform_filter(np.sum(var_inn[:,:,:],axis=0), size=16, mode='constant')
                            #lat[:,:] = lat_rho[:,:]
                            #lon[:,:] = lon_rho[:,:]

                            #value[0,:,:] = 3

                            times[:] = nc.date2num(ens_date, units='days since 1990-01-01',
                             calendar='gregorian')
                            times.units = 'days since 1990-01-01'
                            times.calendar='gregorian'

                            ds.close()

                            
                           # if var equals aicen or vicen
                        
                file_handle.close()
    file_ens.close()

    file_cord_handle.close()

    
def set_number_ensembles(prm_file, ens_count):
    file1 = open(prm_file, 'r') 
    Lines = file1.readlines()
    ii = -1
    for l in Lines:
        ii += 1
        if l[0:7] == 'ENSSIZE':
            #print('Fant den!')
            Lines[ii] = 'ENSSIZE = '+str(ens_count)+'\n'
            #l = 'ENSSIZE = '+str(ens_count)+'\n'

    file1 = open('/Users/sindrefritzner/enkf-c/Assimilation/enkf2.prm', "w")
    file1.writelines(Lines)
    file1.close()


def check_dfs_srf(ens_count, diag_file):
    # Check the DFS and SRF values to see if any exceed the recomended values
    # for example if more than 5% of the data has values larger than the limit
    
    file_handle = nc.Dataset(diag_file, mode='r')
    psrf = file_handle.variables['psrf']
    pdfs = file_handle.variables['pdfs']


    srf_lim = 2 #'må sjekke'
    dfs_lim = ens_count/3 #'må sjekke'
    
    update_srf=np.zeros(psrf.shape[0])
    update_dfs=np.zeros(psrf.shape[0])
    for i in range(0, psrf.shape[0]):
        # more than 5% of the data larger than the limit then reduce
        if len(psrf[i,:,:][psrf[i,:,:]>srf_lim])/(psrf.shape[1]*psrf.shape[2]) > 0.05:
            print('Must increase Rfactor for obstype ' + str(i))
            update_srf[i] = 1
        elif len(psrf[i,:,:][psrf[i,:,:]<srf_lim*0.25])/(psrf.shape[1]*psrf.shape[2]) > 0.99:
            print('Must decrease R_factor for obstype ' + str(i))
            update_srf[i] = -1

        if len(pdfs[i,:,:][pdfs[i,:,:]>dfs_lim])/(pdfs.shape[1]*pdfs.shape[2]) > 0.05:
            print('Must reduce loc_rad for obstype ' + str(i))
            update_dfs[i] = 1
        elif len(pdfs[i,:,:][pdfs[i,:,:]<srf_lim*0.25])/(pdfs.shape[1]*pdfs.shape[2]) > 0.99:
            print('Can increase loc_rad for obstype ' + str(i))
            update_dfs[i] = -1

    # If any must be updated run update prm with information regarding which obstypes that should be updated
    # and how they should be updated. Consider to change back to default values when the assimilation is finished 
    # with satisfying results
    file_handle.close()
    
    return update_srf, update_dfs
    
    
def update_tuning(tuning_file, update_srf, update_dfs):
# To check each locrad and R_factor obstypes.prm should be searched from top to bottom and each time it passes a name
# this name should be rembered such that the next locrad and rfactor encoutered belongs to this name
# It is also important that each obstypes has a locrad specified and an rfractor specified.
# Must investigate a bit more how this should be done in practice, have a list with 0 and 1 per obs type is probably
# the easiest.

    file1 = open(tuning_file, 'r') 
    Lines = file1.readlines()
    ii = -1
    obs_num = -1

    for l in Lines:
        ii += 1
        
        if l[0:4] == 'NAME':
            print(l[7:])
            current_obs = l[7:]
            obs_num += 1
            print(obs_num)
            if obs_num < len(update_srf):
                if update_dfs[obs_num] == 1 and l[0:7] == 'RFACTOR':
                    rf_old = float(l[10:-1])
                    Lines[ii] = 'RFACTOR = '+str(round(rf_old*1.5))+'\n'
                elif update_dfs[obs_num] == -1 and l[0:7] == 'RFACTOR':
                    rf_old = float(l[10:-1])
                    Lines[ii] = 'RFACTOR = '+str(round(rf_old*0.75))+'\n'
                elif update_dfs[obs_num] == -1 and l[0:6] == 'LOCRAD':
                    lr_old = float(l[9:-1])
                    Lines[ii] = 'RFACTOR = '+str(round(lr_old*0.75))+'\n'
                elif update_dfs[obs_num] == -1 and l[0:6] == 'LOCRAD':
                    lr_old = float(l[9:-1])
                    Lines[ii] = 'RFACTOR = '+str(round(lr_old*0.75))+'\n'

def update_the_ensemble(enkf_c_dir, EnKF_var,ens_out_dir,ens_date):
    # Update the ensemble, in practise the whole ensemble does not need to be updated, only those that require 
    # new initial states, but for now everything can be updated. values should also be checked for consistence, 
    # potential large errors should possibly lead to an error, or at least they should be flagged for investigation.

    # Get the list of files that was used as input in the correct order, this file should ideally be written to disk.
    file_ens = open(enkf_c_dir+'files_in_ensemble', 'r') 
    Lines = file_ens.readlines()
    file_count = 0
    prescripts = ['iced.','ocean.']
    syear = str(ens_date.year)
    smnd = str(ens_date.month)
    if ens_date.month < 10:
        smnd = '0' + smnd
    sday = str(ens_date.day)
    if ens_date.day < 10:
        sday = '0' + sday

    zero_checks_0 = ['alvl','qice001','qice002','qice003','qice004','qice005','qice006','qice007'
                    'qsno001','sice001','sice002','sice003','sice004','sice005','sice006','sice007',
                    'vlvl','vsnon']
    qices = ['qice001','qice002','qice003','qice004','qice005','qice006','qice007']
    sices = ['sice001','sice002','sice003','sice004','sice005','sice006','sice007']

    # Make sure that aicen is first in the list of variables such that this can be used for the other variables  
    if EnKF_var.index("aicen") != 0:
        # Switch aicen with the variable that is first in the list
        tlist = copy.deepcopy(EnKF_var)
        tlist[0] = EnKF_var[EnKF_var.index("aicen")]
        tlist[EnKF_var.index("aicen")] = EnKF_var[0]
        EnKF_var = tlist




    for ll in Lines:
        file_count += 1
        if file_count > 1:
            sys.exsit()
        for pre in prescripts:
            file = ens_out_dir+pre+syear+smnd+sday+'_'+ll[0:-1]+'.nc'        
            print(file)
            org_ds = nc.Dataset(file, 'r+', format='NETCDF4')    
            num = str(file_count)
            halo_cells = False
            if file_count < 10:
                num = '0'+num
            
            for var in org_ds.variables.keys():
                if var in EnKF_var:
                    fn = enkf_c_dir+'ensemble_6565/mem0'+num+'_'+var+'.nc.analysis'
                    mem_ds = nc.Dataset(fn, 'r', format='NETCDF4')
                    #print(fn)
                    #print(file)
                    new_var = mem_ds.variables[var]
                    old_var = org_ds.variables[var]

                    # Check bounds for this file, but what should the bounds be? 
                    # With the dfs and srf checks it is not expected that the updates are too large, but could probalby
                    # check just to make sure.
                    # At least SST should never be less than -2 and ice conc should be between 0 and 1.

                    
                    #print(new_var.shape)
                    #print(old_var.shape)

                    if old_var.shape[2]>new_var.shape[3] and pre == 'iced.':
                        halo_cells = True
                    temp = new_var[:]
                    if halo_cells:
                        if len(temp.shape) == 3:
                            temp2 = np.zeros((temp.shape[0],temp.shape[1]+2, 
                                    temp.shape[2]+2))
                            temp2[0,1:-1,1:-1] = temp
                            temp = temp2
                        elif len(temp.shape) == 4:
                            temp2 = np.zeros((temp.shape[0],temp.shape[1], 
                                    temp.shape[2]+2, temp.shape[3]+2))
                            temp2[0,:,1:-1,1:-1] = temp
                            temp = temp2

                    if var == 'temp':
                        # Temperature cannot be below minus 2    
                        temp[temp < -2] = -2
                    if var == 'salt':
                        # Salinity cannot be less than 0    
                        temp[temp < 0] = 0

                    if var == 'aicen':                        
                        # Check ice boundaries
                        temp[temp < 0] = 0
                        temp[temp > 1] = 1
                        # Check that aggregated concentraion is less than 1
                        temsum = np.sum(temp,axis=1)
                        temsum[temsum < 1] = 1
                        
                        for i in range(temp.shape[1]):
                            #print(temsum.shape)
                            #print(temp[:,i,:,:].shape)
                            temp[:,i,:,:] = temp[:,i,:,:]/temsum
                        #print(temp[0,:,:,:].shape)
                        #print(old_var[:].shape)

                        #old_var[:] = temp[0,:,:,:]

                        # These varaibles are to be used for checking the other ice variables,
                        # especially zero checks and new ice checks
                        aicen = temp[:,:,:,:]
                        aice = np.sum(temp,axis=1)
                        aice1 = np.ceil(aice)
    
                    # Set all variables in zero_checks to zero if there is no ice.
                    if var in zero_checks_0:
                        temp[aicen == 0] = 0

                    if var in qices:
                        # set value of new data to -1e8
                        temp[temp > 0] = 0
                        for i in range(temp.shape[1]):
                            temp[:,i,:,:] = np.minimum(temp[:,i,:,:],aice1*-1e8)

                    if var in sices:
                        temp[temp < 0] = 0
                        # set value of new data to 4
                        for i in range(temp.shape[1]):
                            temp[:,i,:,:] = np.maximum(temp[:,i,:,:],aice1*4)
                        

                    if var == 'vicen':                        
                        # Check thickness boundaries
                        temp[temp < 0] = 0
                        # Set thickness to zero for areas without ice
                        temp[aicen == 0] = 0
            
                        # Make sure that new ice is also updated if missed by the assimilation
                        # Assume that the new thickness is very thin, vicen=aicen just for simplicity,
                        # this is not really expected to happen, but can cause numerical errors
                        temp = np.maximum(temp,aicen)                        
                       
                    if pre == 'iced.':
                        old_var[:] = temp[0]
                    elif pre == 'ocean.':
                        old_var[:] = temp

                    
                    mem_ds.close()
            org_ds.close()


            #old_var = new_var
        
