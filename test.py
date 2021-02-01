import python_tools.enkf_c_toolbox_parallel as et
import python_tools.Obs_sindre as ob
from datetime import datetime
import netCDF4 as nc
import numpy as np
import sys
from shutil import copyfile
import os

Prep = False
Assimilate = False
Update = False
Write = False
Obs = False
ens_count = False

mode = input("What to do?[Prep,Assimilate,Update,Write,Obs,All]: ") 
if mode == 'Prep':
    Prep = True
elif mode == 'Assimilate':
    Assimilate = True
elif mode == 'Update':
    Update = True
elif mode == 'Write':
    Write = True
elif mode == 'Obs':
    Obs = True
elif mode == 'ens_count':
    ens_count = True
elif mode == 'All':
    Prep = True
    Assimilate = True
    Update = True 
    Write = True
    Obs = True
else:
    sys.exit(mode+' is not implemented yet')

computer = input("Type in computer in use[met_local(default),nebula,fram]: ") or "met_local" 
print(computer)
end_date = datetime(2018,1,8)
if computer == 'met_local':
    res_dir = "/home/sindremf/PHD2/Work/Test_assimiation/Resdir/"
    grid_dir = '/home/sindremf/PHD2/Work/Barents/data_dir/org_files/barents_grd.nc'
    enkf_c_dir = "/home/sindremf/PHD2/Work/Assim_enkf-c/"
    obs_dir = '/home/sindremf/PHD2/Work/Observations/'
    Nens = 10 # This is maximum value, but not neccesarily the number used.
elif computer == 'nebula':
    res_dir = '/nobackup/forsk/sm_sinfr/Results/barents/'
    grid_dir = '/home/sm_sinfr/metroms_apps/barents-2.5km/grid/barents_grd.nc'
    enkf_c_dir = "/nobackup/forsk/sm_sinfr/Assim_enkf-c/"
    obs_dir = '/nobackup/forsk/sm_sinfr/Observations/'
    save_dir = '/nobackup/forsk/sm_sinfr/Results/Assim_res/'
    Nens = 10 # This is maximum value, but not neccesarily the number used.
elif computer == 'fram':
    res_dir = '/cluster/work/users/sfr009/Results/barents/'
    grid_dir = '/cluster/home/sfr009/metroms_apps/barents-2.5km/grid/barents_grd.nc'
    enkf_c_dir = '/cluster/work/users/sfr009/Assim_enkf-c/'
    obs_dir = '/cluster/work/users/sfr009/Observations/'
    save_dir = '/cluster/work/users/sfr009/Results/Assim_res/'
    Nens = 10 # This is maximum value, but not neccesarily the number used.
else:
    sys.exit(computer+' not yet implementet')

res_type = 'ice'
EnKF_var = ['aicen','vicen','temp','salt','qice001','qice002','qice003','qice004',
            'qice005','qice006','qice007','qsno001','sice001','sice002','sice003','sice004',
            'sice005','sice006','sice007','vsnon', 'Tsfcn']

# Write which and how many ensemble members that are ready

# this one is only for testing, more correct to do this within prep in case more ensembles are
# finishing during assimilation
if ens_count:
    ens_count = et.count_ens(date=end_date,enkf_c_dir=enkf_c_dir,res_dir=res_dir)

# Prep the ensemble
if Prep:
    #Return Nens!!!
    et.Prep_ensemble(ens_date = end_date, grid_dir=grid_dir, 
                     ens_inn_dir=res_dir, enkf_c_dir =enkf_c_dir, 
                     res_type = 'ice', EnKF_var=EnKF_var,
                     Nens = Nens)
    ens_count = et.count_ens(date=end_date,enkf_c_dir=enkf_c_dir,res_dir=res_dir)

# copy and prep observations
# See also Keguang code for download
if Obs:
    ob.prep_osisaf_obs(date=end_date,obs_dir=obs_dir+'/Org_OSISAF/', Assim_dir=enkf_c_dir)
    

    # Copy the AMSR observation from obs directory
    file_amsr_inn = obs_dir+'AMSR_observations/'+'amsr2_'+end_date.strftime('%Y%m%d')+'.nc'
    file_amsr_out = enkf_c_dir+'obs/AMSR/this_day.nc'
    if os.path.exists(file_amsr_inn):
        copyfile(file_amsr_inn,file_amsr_out)
    else:
        print('AMSR file not availible: '+file_amsr_inn)
        
        try:
            os.remove(file_amsr_out)
        except:
            pass
        
    

# run the assimlation
if Assimilate:
    copyfile(enkf_c_dir+'bld/'+'Makefile.'+computer,enkf_c_dir+'Makefile')
    et.cmd('make enkf')

#Update the ensemble
if Update:
    et.update_the_ensemble(enkf_c_dir =enkf_c_dir, EnKF_var=EnKF_var, ens_out_dir=res_dir,ens_date = end_date)


#Write important results to file.
if Write:
    # Kun for testing! Blir lagret etter prep i det operasjonelle!
    file_ens = open(enkf_c_dir+'files_in_ensemble', 'r') 
    Lines = file_ens.readlines()
    ens_count = len(Lines)
    ###############################
    et.write_results(date=end_date,enkf_c_dir=enkf_c_dir,ens_out_dir=res_dir, Nens=ens_count, save_dir=save_dir)
