clear all

kmt = ncread('kmt.nc','kmt');
lon = ncread('grid.nc','ulon');
lat = ncread('grid.nc','ulat');

depth_t = ncread('A20_grd_openBering.nc','h');

[a,b] = size(lon);

c = 5;

for i = 1:c
    zt(i) = (1 + 2*i)*2*i;
end

num_levels = c*ones(a,b);

num_levels(kmt == 0) = 0;

file = 'grid_spec.nc';
delete(file);

lon = lon + 180;
nccreate(file,'lon','Dimensions',{'i',a,'j',b})
ncwrite(file,'lon',lon) 

nccreate(file,'lat','Dimensions',{'i',a,'j',b})
ncwrite(file,'lat',lat) 

nccreate(file,'zt','Dimensions',{'k',c})
ncwrite(file,'zt',zt) 

nccreate(file,'depth_t','Dimensions',{'i',a,'j',b})
ncwrite(file,'depth_t',depth_t) 

nccreate(file,'num_levels','Dimensions',{'i',a,'j',b})
ncwrite(file,'num_levels',num_levels) 
