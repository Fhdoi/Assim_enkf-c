# Generate forcingfiles requires from MET Barents
import os
import glob
import shutil
import argparse
import subprocess
import re
from copy import copy
from datetime import datetime, timedelta
from subprocess import call
import netCDF4 as nc
import numpy as np
#import wrf # Denne har compabilitets problemer
import sys

def cmd(command,ii=1, exits=True):
    """Function runs provided command in the system shell.

    Args:
        command (string) : Command to be executed in shell
    Returns:
        result (integer) : Shell returned status of command
    """
    print("> " + command)
    result = subprocess.call(command, shell = True)

    if result != 0:
        
        print("Command failed, "+str(ii)+"time: %d" % result)
        ii += 1
        if ii > 5:
            sys.exit('This command dosen\'t work: '+command)
        if exits:
            cmd(command,ii)
        else:
            return result
        
        
    else:
        return result
        
def replace_line(filename, pattern, replacer, precursor="", replace_all=False):
    """
    Function that searches a text file line-by-line and replaces a line with user speified string
    if a user provided regex pattern is matched on that line. The search only starts after an
    optional precursor string is found in the file.

    Args:
        filename (str)     : Filename of text file
        pattern (str)      : Regular expression pattern to search for
        replacer (str)     : String to replace on line where regex was matched
        precursor (str)    : Optional string to limit the regex search start from
                             the line where the precursor was found. If not provided,
                             search starts from the top of file.
        replace_all (bool) : If True, replace all lines where regex matched (after
                             optional precursor), otherwise replace only first occurence
    """
    tmp_file = filename + "_tmp"
    begin_search = False
    done_replace = False

    with open(tmp_file, "w") as new_file:
        with open(filename, "r") as old_file:
            for line in old_file:
                if precursor in line:
                    begin_search=True  # only look for line after precursor

                # if precursor found, no
                if begin_search and not done_replace and re.match(pattern, line):
                    new_file.write(replacer)

                    if not replace_all:
                        done_replace = True
                else:
                    new_file.write(line)  # otherwise just copy existing line

    shutil.move(tmp_file, filename)  # replace original with newly written
    
def valid_date(date_string):
    """
    Function that validates if a date string is of a certain format
    and returns a datetime object corresponding to that string.

    Args:
        date_string (str) : String of a date to be checked
    Returns:
        datetime (datetime) : Datetime object created from input string
    """
    try:
        return datetime.strptime(date_string, '%Y%m%d')
    except ValueError:
        message = 'Not a valid date: {0}.'.format(date_string)
        raise argparse.ArgumentTypeError(message)
        
def valid_path(path):
    """
    Function that validates if path exists.

    Args:
        path (str) : Suggested file path
    Returns:
        path (str) : If path exists
    """
    if os.path.exists(path):
        return path
    else:
        raise argparse.ArgumentTypeError("Invalid path {}!".format(path))
        
def fetch_clm_HI(cfg_file, tmp_topaz2,start_date,work_dir,numdays=10):
    # Delete the temp file if it exists
    end_date = start_date + timedelta(days=numdays)
    try:
        os.remove(tmp_topaz2)
    except:
        pass
    # replace whatever output file is in cfg file with a tmp filename
    replace_line(cfg_file, "file ?= ?.+", "file={}\n".format(tmp_topaz2), precursor="[output]")
    print('\nFetching TOPAZ data from {} to {} and storing data in {}...\n'.format(start_date, end_date, work_dir))
    #cmd('/usr/bin/fimex-1.5 -c {} --extract.reduceTime.start {} --extract.reduceTime.end {} --input.config {}'.format(cfg_file, start_date.strftime('%Y-%m-%d'), end_date.strftime('%Y-%m-%d'), ncml_file))
    cmd('/usr/bin/fimex-1.5 -c {} --extract.reduceTime.start {} --extract.reduceTime.end {}'.format(cfg_file, start_date.strftime('%Y-%m-%d'), end_date.strftime('%Y-%m-%d')))
    if os.path.exists(tmp_topaz2):
            print('Successfull horizontal interpolation!')
    else:
            raise RuntimeError('Failed horizontal interpolation! Aborting...')
    tpz_nlevs = nc.Dataset(tmp_topaz2).dimensions['depth'].size
    tmpnc_filepath = os.path.join(work_dir, 'tmp.nc')
    cmd('cp {} {}'.format(tmp_topaz2, tmpnc_filepath))
    fid = nc.Dataset(tmp_topaz2)
    for key in fid.variables.keys():
        if '_FillValue' in fid.variables[key].ncattrs():
            cmd('ncatted -O -h -a _FillValue,{},m,f,0.0 {}'.format(key, tmpnc_filepath))
        else:
            print('Skipping removal of _FillValue for '+key)
    fid.close()
    
    cmd('mv {} {}'.format(tmpnc_filepath, tmp_topaz2))
    
def generate_empty_clmfile(clmfile='barents_clm.nc',roms_nlevs=42,roms_nx=949,roms_ny=739):
    clim_var2d = ['zeta','ubar','vbar','aice','hice','uice','vice','snow_thick','ageice']
    clim_var3d = ['u','v','salt','temp']

    try:
        os.remove(clmfile)
    except:
        pass
    ds = nc.Dataset(clmfile, 'w', format='NETCDF4')


    time = ds.createDimension('clim_time')

    # Jeg er veldig usikker paa hvorfor det er slik, bor sporre Jostein/Johannes
    z = ds.createDimension('s_rho', roms_nlevs)
    x = ds.createDimension('eta_u', roms_nx)
    y = ds.createDimension('xi_u', roms_ny-1)

    x = ds.createDimension('eta_v', roms_nx-1)
    y = ds.createDimension('xi_v', roms_ny)

    x = ds.createDimension('eta_rho', roms_nx)
    y = ds.createDimension('xi_rho', roms_ny)



    times = ds.createVariable('clim_time', 'f4', ('clim_time'))

    for i in clim_var2d:
        if i == 'ubar':
            ds.createVariable(i, 'f4', ('clim_time','eta_u','xi_u',))
        elif i == 'vbar':
            ds.createVariable(i, 'f4', ('clim_time','eta_v','xi_v',))
        else:
            ds.createVariable(i, 'f4', ('clim_time','eta_rho','xi_rho',))
        
    for i in clim_var3d:
        if i == 'u':
            ds.createVariable(i, 'f4', ('clim_time','s_rho','eta_u','xi_u',))
        elif i == 'v':
            ds.createVariable(i, 'f4', ('clim_time','s_rho','eta_v','xi_v',))
        else:
            ds.createVariable(i, 'f4', ('clim_time','s_rho','eta_rho','xi_rho',))
        
                
        
    ds.close()


def fill_clim(tmp_topaz2, clmfile='barents_clm.nc', grid_file='barents_grd.nc', roms_nx=949,roms_ny=739, roms_nlevs=42):
    sigma = np.flip(np.arange(1/(roms_nlevs*2),1,1/roms_nlevs))
    file_grid = nc.Dataset(grid_file)
    geo_depth = file_grid['h'][:].data
    file_grid.close()
    tpz_nlevs = nc.Dataset(tmp_topaz2).dimensions['depth'].size
    nots = ['x','y','depth','time','model_depth','projection_lcc','lon','lat']
    clim_var2d = ['zeta','aice','hice']
    tpz = nc.Dataset(tmp_topaz2)
    ds = nc.Dataset(clmfile, 'r+')
    ctime = ds['clim_time']
    #ctime[:] = tpz['time'][:].data
    #ctime.units = tpz['time'].units
    ctime[:] = nc.date2num(nc.num2date(tpz['time'][:].data,tpz['time'].units),'days since 1970-1-1T00:00:00Z')
    ctime.units = 'days since 1970-1-1T00:00:00Z'

    tpz_depth=tpz['depth'][:].data
    tpz_sigma = np.zeros((tpz_nlevs,roms_nx,roms_ny))
    temp = np.transpose(np.tile(tpz_depth.data,[roms_ny,roms_nx,1]))
    #depth2[depth2 < 100.0] = 3000
    for i in range(0,12):
        tpz_sigma[i,:,:] = temp[i,:,:]/geo_depth
        

        
    # Zeropad tpz_sigma
    tpz_sigma = np.append(np.zeros((1,roms_nx,roms_ny)),tpz_sigma,axis=0)

    # Fix problem when ocean is deeper than maximum tpz depth
    bottom_sigma = tpz_sigma[-2:-1,:,:]
    # make sure that there are values below to interpolate, does not effect the result or provide any more information,
    # just helps numerically to avoid nan values for the interpolation
    bottom_sigma[bottom_sigma < 1] = 1 
    tpz_sigma = np.append(tpz_sigma,bottom_sigma,axis=0)



    for i in tpz.variables.keys():
    #for i in ['u']:
        if i not in nots:
            temp = ds[i]
            # Sjekk if we are writing a 2 or 3 dimensional variable
            if len(temp.shape) == 3:
                if i == 'ubar':
                    temp[:] = tpz[i][:,:,:-1].data
                elif i == 'vbar':
                    temp[:] = tpz[i][:,:-1,:].data
                else:
                    temp[:] = tpz[i][:].data
                
            elif len(temp.shape) == 4:
                # Do vertical interpolation
                test = tpz[i]
                # surface and bottom pad (This helps the interpolation)
                test = np.append(np.append(test[:,0:1,:,:],test,axis=1),test[:,-2:-1,:,:],axis=1)
                # Loop through the dates in the file
                #print(np.sum(np.isnan(test)))
                #print(np.sum(test > 1000))
                test3d = np.zeros((test.shape[0],roms_nlevs, roms_nx, roms_ny))
                for j in range(0,temp.shape[0]):            
                    # add surface value for interpolation
                    #print(np.sum(np.isnan(temp)))
                    test3d[j,:,:,:] = wrf.interplevel(test[j,:,:,:], tpz_sigma, sigma, missing=0.0, meta=False)
                
                # I think missing should fix this, but it does not seem to work
                #test3d = temp2[:,:,:,:]
                test3d[test3d > 1000] = 0
                if i == 'u':
                    temp[:] = test3d[:,:,:,:-1]
                elif i == 'v':
                    temp[:] = test3d[:,:,:-1,:]
                else:
                    temp[:] = test3d 
                
            temp.time = 'clim_time'
                
                
                # if i equal u or v also make ubar
                
            #print('etter')
            #print(np.sum(np.isnan(temp)))
            #print(np.sum(temp[:].data > 1000))
    

    # Add ubar and vbar to the clm file
    #ds = nc.Dataset(clmfile, 'r+')
    u = ds['u']
    v = ds['v']
    ubar = ds['ubar']
    vbar = ds['vbar']
    dsigma = 1/roms_nlevs

    ubar[:] = np.sum(dsigma*u[:],axis=1)
    vbar[:] = np.sum(dsigma*v[:],axis=1)

    ubar.time = 'clim_time'
    vbar.time = 'clim_time'
    
    
            
    ds.close()
    tpz.close()
    
    
def generate_bry(clmfile='barents_clm.nc', bryfile='barents_bry.nc',roms_nlevs=42, roms_nx=949, roms_ny=739):
    # Generate bry file

    cl = nc.Dataset(clmfile, 'r')

    try:
        os.remove(bryfile)
    except:
        pass
    ds = nc.Dataset(bryfile, 'w', format='NETCDF4')


    time = ds.createDimension('bry_time')
    z = ds.createDimension('s_rho', roms_nlevs)
    x = ds.createDimension('eta_rho', roms_nx)
    y = ds.createDimension('xi_rho', roms_ny)
    y2 = ds.createDimension('xi_u', roms_ny-1)
    x2 = ds.createDimension('eta_v', roms_nx-1)



    times = ds.createVariable('bry_time', 'f4', ('bry_time'))
    times[:] = cl['clim_time'][:].data
    times.units = cl['clim_time'].units


    for i in cl.variables.keys():
        if i not in ['clim_time']:
            clm_var = cl[i]
            if len(clm_var.shape) == 3:
                if i == 'ubar':
                    temp_north = ds.createVariable(i+'_north', 'f4', ('bry_time','xi_u',))
                    temp_south = ds.createVariable(i+'_south', 'f4', ('bry_time','xi_u',))
                    temp_west = ds.createVariable(i+'_west', 'f4', ('bry_time','eta_rho',))
                    temp_east = ds.createVariable(i+'_east', 'f4', ('bry_time','eta_rho',))

                    temp_north[:] = clm_var[:,-1,:]
                    temp_south[:] = clm_var[:,0,:]
                    temp_west[:] = clm_var[:,:,0]
                    temp_east[:] = clm_var[:,:,-1]
                    
                elif i == 'vbar':
                    temp_north = ds.createVariable(i+'_north', 'f4', ('bry_time','xi_rho',))
                    temp_south = ds.createVariable(i+'_south', 'f4', ('bry_time','xi_rho',))
                    temp_west = ds.createVariable(i+'_west', 'f4', ('bry_time','eta_v',))
                    temp_east = ds.createVariable(i+'_east', 'f4', ('bry_time','eta_v',))

                    temp_north[:] = clm_var[:,-1,:]
                    temp_south[:] = clm_var[:,0,:]
                    temp_west[:] = clm_var[:,:,0]
                    temp_east[:] = clm_var[:,:,-1]
                    
                else:
                    temp_north = ds.createVariable(i+'_north', 'f4', ('bry_time','xi_rho',))
                    temp_south = ds.createVariable(i+'_south', 'f4', ('bry_time','xi_rho',))
                    temp_west = ds.createVariable(i+'_west', 'f4', ('bry_time','eta_rho',))
                    temp_east = ds.createVariable(i+'_east', 'f4', ('bry_time','eta_rho',))

                    temp_north[:] = clm_var[:,-1,:]
                    temp_south[:] = clm_var[:,0,:]
                    temp_west[:] = clm_var[:,:,0]
                    temp_east[:] = clm_var[:,:,-1]
            
            elif len(clm_var.shape) == 4:
                if i == 'u':
                    temp_north = ds.createVariable(i+'_north', 'f4', ('bry_time','s_rho','xi_u',))
                    temp_south = ds.createVariable(i+'_south', 'f4', ('bry_time','s_rho','xi_u',))
                    temp_west = ds.createVariable(i+'_west', 'f4', ('bry_time','s_rho','eta_rho',))
                    temp_east = ds.createVariable(i+'_east', 'f4', ('bry_time','s_rho','eta_rho',))

                    temp_north[:] = clm_var[:,:,-1,:]
                    temp_south[:] = clm_var[:,:,0,:]
                    temp_west[:] = clm_var[:,:,:,0]
                    temp_east[:] = clm_var[:,:,:,-1]
                    
                elif i == 'v':
                    temp_north = ds.createVariable(i+'_north', 'f4', ('bry_time','s_rho','xi_rho',))
                    temp_south = ds.createVariable(i+'_south', 'f4', ('bry_time','s_rho','xi_rho',))
                    temp_west = ds.createVariable(i+'_west', 'f4', ('bry_time','s_rho','eta_v',))
                    temp_east = ds.createVariable(i+'_east', 'f4', ('bry_time','s_rho','eta_v',))

                    temp_north[:] = clm_var[:,:,-1,:]
                    temp_south[:] = clm_var[:,:,0,:]
                    temp_west[:] = clm_var[:,:,:,0]
                    temp_east[:] = clm_var[:,:,:,-1]
                    
                else:
                    temp_north = ds.createVariable(i+'_north', 'f4', ('bry_time','s_rho','xi_rho',))
                    temp_south = ds.createVariable(i+'_south', 'f4', ('bry_time','s_rho','xi_rho',))
                    temp_west = ds.createVariable(i+'_west', 'f4', ('bry_time','s_rho','eta_rho',))
                    temp_east = ds.createVariable(i+'_east', 'f4', ('bry_time','s_rho','eta_rho',))

                    temp_north[:] = clm_var[:,:,-1,:]
                    temp_south[:] = clm_var[:,:,0,:]
                    temp_west[:] = clm_var[:,:,:,0]
                    temp_east[:] = clm_var[:,:,:,-1]
            
            temp_north.time = 'bry_time'
            temp_south.time = 'bry_time'
            temp_west.time = 'bry_time'
            temp_east.time = 'bry_time'

        
    cl.close()
        
    ds.close()
    
def download_arome(start_date,work_dir,final_atm_file,cfg_file='fimex_barents-2.5km_arome.cfg',numdays=10,hourly=False):
    odir = 'https://thredds.met.no/thredds/dodsC/aromearcticarchive/'
    otype = '/arome_arctic_extracted_2_5km_'
    #final_atm_file = work_dir + 'barents_atm.nc'
    date_list = [start_date + timedelta(days=x) for x in range(numdays)]
    filelist = []
    for date in date_list:
        syear = str(date.year)
        smnd = str(date.month)
        sday = str(date.day)
        if date.month < 10:
            smnd = '0'+smnd
        if date.day < 10:
            sday = '0'+sday
        atmfile=work_dir + 'barents_atm_'+syear+smnd+sday+'.nc'
        filelist.append(atmfile)
        # replace whatever output file is in cfg file with a tmp filename
        otimes = ['00','06','12','18']
        for i in range(0,4):
            otime = otimes[i]
            input_file = odir+'/'+syear+'/'+smnd+'/'+sday+otype+syear+smnd+sday+'T'+otime+'Z.nc'
            output_file = work_dir+'tmp.nc'
            replace_line(cfg_file, "file ?= ?.+", "file={}\n".format(input_file), precursor="[input]")
            replace_line(cfg_file, "file ?= ?.+", "file={}\n".format(output_file), precursor="[output]")
            
            Result =cmd('/usr/bin/fimex-1.5 -c '+cfg_file,exits=False)
            file_last = False

            # If data did not exist, use previous time step
            if Result != 0:
                file_last = True
                num = previous_time(date=date,i=i,odir=odir,otype=otype,work_dir=work_dir,cfg_file=cfg_file)
                #print('Dette skjer!')
                #if i == 0:
                #    otime2 = otimes[3]
                #    # denne fungerer foreløpig kun for dager større enn 10
                    # Og bør også være iterativ hvis det mangler flere på rad
                    # Muligens letter å bare bruke den Jostein bruker
                #    sday = str(date.day-1)    
                #else:
                #    otime2 = otimes[i -1]
                #input_file = odir+'/'+syear+'/'+smnd+'/'+sday+otype+syear+smnd+sday+'T'+otime2+'Z.nc'
                #output_file = work_dir+'tmp.nc'
                #replace_line(cfg_file, "file ?= ?.+", "file={}\n".format(input_file), precursor="[input]")
                #replace_line(cfg_file, "file ?= ?.+", "file={}\n".format(output_file), precursor="[output]")
                #file_last = True
                #cmd('/usr/bin/fimex-1.5 -c '+cfg_file,exits=True)
        

            if hourly:
                #Every hour
                cmd('ncks -O -d time,0,5 '+work_dir+'tmp.nc '+work_dir+'tmp'+otime+'.nc')
            else:
                #Every 6 hours
                cmd('ncks -O -d time,0 '+work_dir+'tmp.nc '+work_dir+'tmp'+otime+'.nc')
            # If previous day was used due to missing data
            if file_last:
                cmd('ncap2 -O -h -s "time=time + (3600*6*'+str(num)+')" '+work_dir+'tmp'+otime+'.nc'+' '+work_dir+'tmp'+otime+'.nc')
            
            #
        cmd('ncrcat -O '+work_dir+'tmp00.nc '+work_dir+'tmp06.nc '+work_dir+'tmp12.nc '+work_dir+'tmp18.nc '+atmfile)

        # Remove unnecessary dimensions
        # Denne fungerer kun i commando vindu av en eller annen grunn:/
        # Conda Enviroment problem!!!, fungerer ikke med Barents lenger, men bruk bare test
        cmd('ncwa -O -h -a height1,height3,height_above_msl,height0 '+atmfile+' '+ atmfile)
        
        cmd('ncatted -O -h -a units,Pair,m,c,millibar '+atmfile)
        cmd('ncap2 -O -h -s "Pair=0.01*Pair" '+atmfile+' '+atmfile)

        # Fix time 
        cmd('ncap2 -O -h -s "time=time/(3600*24.0)" '+atmfile+' '+atmfile)

        # Fix Tair
        cmd('ncatted -O -h -a units,Tair,m,c,Celsius '+atmfile)
        cmd('ncap2 -O -h -s "Tair=Tair-273.16" '+atmfile+' '+atmfile)

        # Fix Qair
        cmd('ncatted -O -h -a units,Qair,m,c,percentage '+atmfile)
        cmd('ncap2 -O -h -s "Qair_rel=100.0*Qair_rel" '+atmfile+' '+atmfile)

        # Fix rain, data is precip./hour
        cmd('ncatted -O -h -a units,rain,m,c,"kilogram meter-2 second-1" '+atmfile)
        cmd('ncap2 -O -h -s "rain=rain/3600.0" '+atmfile+' '+atmfile)

        # Fix Uwind/Vwind attributes
        cmd('ncatted -O -h -a units,Uwind,m,c,"meter second-1" '+atmfile)
        cmd('ncatted -O -h -a units,Vwind,m,c,"meter second-1" '+atmfile)

        
    cmd('ncrcat -O '+convert_list_to_string(filelist)+' '+final_atm_file)
    
    # Dette kan gjøres når hele perioden er lastet ned
    # Minne problem, gjøre det for hver dag tror jeg
    # Fix Pair
    #cmd('ncatted -O -h -a units,Pair,m,c,millibar '+final_atm_file)
    #cmd('ncap2 -O -h -s "Pair=0.01*Pair" '+final_atm_file+' '+final_atm_file)

    # Fix time 
    #cmd('ncap2 -O -h -s "time=time/(3600*24.0)" '+final_atm_file+' '+final_atm_file)

    # Fix Tair
    #cmd('ncatted -O -h -a units,Tair,m,c,Celsius '+final_atm_file)
    #cmd('ncap2 -O -h -s "Tair=Tair-273.16" '+final_atm_file+' '+final_atm_file)

    # Fix Qair
    #cmd('ncatted -O -h -a units,Qair,m,c,percentage '+final_atm_file)
    #cmd('ncap2 -O -h -s "Qair=100.0*Qair" '+final_atm_file+' '+final_atm_file)

    # Fix rain, data is precip./hour
    #cmd('ncatted -O -h -a units,rain,m,c,"kilogram meter-2 second-1" '+final_atm_file)
    #cmd('ncap2 -O -h -s "rain=rain/3600.0" '+final_atm_file+' '+final_atm_file)

    # Fix Uwind/Vwind attributes
    #cmd('ncatted -O -h -a units,Uwind,m,c,"meter second-1" '+final_atm_file)
    #cmd('ncatted -O -h -a units,Vwind,m,c,"meter second-1" '+final_atm_file)

def previous_time(date,i,odir,otype,work_dir,cfg_file,num=0):
    num = num+1
    otimes = ['00','06','12','18']
    if i == 0:
        otime = otimes[3]
        date = date + timedelta(days=-1)
            
        sday = str(date.day) if date.day > 9 else '0'+str(date.day)
        smnd = str(date.month) if date.month > 9 else '0'+str(date.month)
        syear = str(date.year)
    else:
        otime = otimes[i -1]

    input_file = odir+'/'+syear+'/'+smnd+'/'+sday+otype+syear+smnd+sday+'T'+otime+'Z.nc'
    output_file = work_dir+'tmp.nc'
    replace_line(cfg_file, "file ?= ?.+", "file={}\n".format(input_file), precursor="[input]")
    replace_line(cfg_file, "file ?= ?.+", "file={}\n".format(output_file), precursor="[output]")
    file_last = True
    Res = cmd('/usr/bin/fimex-1.5 -c '+cfg_file,exits=False)
    if Res != 0:
        previous_time(date,i,odir,otype,work_dir,cfg_file,num)
    else:
        return num
    
    
def convert_list_to_string(org_list, seperator=' '):
    """ Convert list to string, by joining all item in list with given separator.
        Returns the concatenated string """
    return seperator.join(org_list)


