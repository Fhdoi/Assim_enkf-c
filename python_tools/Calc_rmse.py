import netCDF4 as nc
import datetime
import numpy as np
import matplotlib.pyplot as plt

def calc_rmse(var1,var2,kmt):
    return np.sqrt(sum(sum(kmt[:,:]*(var1[:,:] - var2[:,:])**2))/sum(sum(kmt)))

assim_dir = '/home/sindremf/PHD2/Work/Testing/Check_ensemble/'
start_date = datetime.datetime(2018,1,15)
timestep = 7
num_dates = 3
rmse_inn = np.zeros(num_dates)
rmse_out = np.zeros(num_dates)
for j in range(num_dates):
    date = start_date + datetime.timedelta(days=timestep*j)
    syear = str(date.year)
    smnd = str(date.month) if date.month > 9 else '0'+str(date.month)
    sday = str(date.day) if date.day > 9 else '0'+str(date.day)
    file = assim_dir+'Assim_summary_'+syear+smnd+sday+'.nc'
    print(file)
    handle = nc.Dataset(file)
    aice_inn = handle['aice_inn'][:].data
    aice_out = handle['aice_out'][:].data
    obs = handle['Obs1'][:].data

    # Filter land values
    kmt = 1-np.isnan(obs)*1
    obs[np.isnan(obs)] = 0

    # Calculate the average errror over all ensemble members
    rinn = 0
    rout = 0
    for i in range(aice_inn.shape[1]):
        rinn += calc_rmse(aice_inn[0,i,:,:],obs[0,:,:],kmt[0,:,:])
        rout += calc_rmse(aice_out[0,i,:,:],obs[0,:,:],kmt[0,:,:])
    rmse_inn[j] = rinn/aice_inn.shape[1]
    rmse_out[j] = rout/aice_inn.shape[1]
print(rmse_inn)
print(rmse_out)
