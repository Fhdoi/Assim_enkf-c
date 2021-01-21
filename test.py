import enkf_c_toolbox as et
from datetime import datetime
import netCDF4 as nc
import numpy as np

end_date = datetime(2018,1,15)
computer = 'met_local'
if computer == 'met_local'
    res_dir = "/home/sindremf/PHD2/Work/Test_assimiation/Resdir/"
    grid_dir = '/home/sindremf/PHD2/Work/Barents/data_dir/org_files/barents_grd.nc'
    enkf_c_dir = "/home/sindremf/PHD2/Work/Test_assimiation/"
    Nens = 10
elif computer == 'nebula'
    sys.exit(computer+' not yet implementet')
else:
    sys.exit(computer+' not yet implementet')

res_type = 'ice'
EnKF_var = ['aicen','vicen','temp','salt','qice001','qice002','qice003','qice004',
            'qice005','qice006','qice007','sice001','sice002','sice003','sice004',
            'sice005','sice006','sice007']
Prep = True
Assimilate = True
Update = True
Write_res = True



# Prep the ensemble
if Prep:
    et.Prep_ensemble(ens_date = end_date, grid_dir=grid_dir, 
                     ens_inn_dir=res_dir, enkf_c_dir =enkf_c_dir, 
                     res_type = 'ice', EnKF_var=EnKF_var,
                     Nens = Nens)

# copy in observation
#

# run the assimlation
#if Assimilate:
    #et.cmd('make enkf')

#Update the ensemble
if Update:
    et.update_the_ensemble(enkf_c_dir =enkf_c_dir, EnKF_var=EnKF_var, ens_out_dir=res_dir,ens_date = end_date)


#Write important results to file.
#if Write_res:
