import enkf_c_toolbox as et
from datetime import datetime
import netCDF4 as nc
import numpy as np
import sys

Prep = False
Assimilate = False
Update = False
Write_res = False
Obs = False

mode = input("What to do?[Prep,Assimilate,Update,Write_res,Obs,All]: ") 
if mode == 'Prep':
    Prep = True
elif mode == 'Assimilate':
    Assimilate = True
elif mode == 'Update':
    Update = True
elif mode == 'Write_res':
    Write_res = True
elif mode == 'Obs':
    Obs = True
elif mode == 'All':
    Prep = True
    Assimilate = True
    Update = True 
    Write_res = True
    Obs = True
else:
    sys.exit(mode+' is not implemented yet')

computer = input("Type in computer in use[met_local(default),nebula,fram]: ") or "met_local" 
print(computer)
end_date = datetime(2018,1,15)
if computer == 'met_local':
    res_dir = "/home/sindremf/PHD2/Work/Test_assimiation/Resdir/"
    grid_dir = '/home/sindremf/PHD2/Work/Barents/data_dir/org_files/barents_grd.nc'
    enkf_c_dir = "/home/sindremf/PHD2/Work/Assim_enkf-c/"
    obs_dir = '/home/sindremf/PHD2/Work/Observations/'
    Nens = 10
elif computer == 'nebula':
    sys.exit(computer+' not yet implementet')
elif computer == 'fram':
    sys.exit(computer+' not yet implementet')
else:
    sys.exit(computer+' not yet implementet')

res_type = 'ice'
EnKF_var = ['aicen','vicen','temp','salt','qice001','qice002','qice003','qice004',
            'qice005','qice006','qice007','sice001','sice002','sice003','sice004',
            'sice005','sice006','sice007']




# Prep the ensemble
if Prep:
    et.Prep_ensemble(ens_date = end_date, grid_dir=grid_dir, 
                     ens_inn_dir=res_dir, enkf_c_dir =enkf_c_dir, 
                     res_type = 'ice', EnKF_var=EnKF_var,
                     Nens = Nens)

# copy and prep observations
# See also Keguang code for download
if Obs:
    et.get_osisaf_obs(date=end_date,obs_dir=obs_dir+'/Org_OSISAF/', Assim_dir=enkf_c_dir)

# run the assimlation
if Assimilate:
    et.cmd('make enkf')

#Update the ensemble
if Update:
    et.update_the_ensemble(enkf_c_dir =enkf_c_dir, EnKF_var=EnKF_var, ens_out_dir=res_dir,ens_date = end_date)


#Write important results to file.
#if Write_res:
