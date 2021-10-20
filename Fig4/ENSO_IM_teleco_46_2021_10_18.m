% 03.05.21. This is to get data from netcdf files and save them in a mat
% file, to be used by further "computation scripts". Companion to *_45.m.
% Instead of saving monthly means, i save seasonal means, because monthly
% means would result in large difficult-to-open mat files. Already with 50
% members, i had files of 10GB containing matrices of
% (200*300*12*250*50=9e9) entries per field. Also, i just get Nino34.
% 04.05.21. Luckily, there is now output of regridded SST in the atm
% directories. With this, the code is a lot simplified and much faster too
% (as i don't need to load the 3D var TEMP); the parfor loops took only ~10
% minutes... but saving is another ~20 mins. Anyhow, because of just
% looking at the Nino34 box, i don't need the landmask.

% 18.10.21. Tidy code.

% This code is shared as supporting material for the paper: 
% 
% “Ubiquity of human-induced changes in climate variability”, Rodgers et
% al., Earth System Dynamics, 2021. https://doi.org/10.5194/esd-2021-50
%
% specifically, for Figs. 4 and S9. The order in which the scripts can be
% used is: 
% 1. ENSO_IM_teleco_46_2021_10_18.m 
% 2. ENSO_IM_teleco_45_2021_10_18.m 
% 3. ENSO_IM_teleco_41_2021_10_18.m
% 4. ENSO_IM_teleco_presentation_2021_10_18.m

% Author: Tamas Bodai (bodai at pusan.ac.kr, bodait at yahoo.com)

clear 
close all

% Add libraries to path if not already in there
lib =    {'/home/zmaw/u234078/Documents/MATLAB/Useful'; ...
                '/Users/bodai/Documents/MATLAB/Useful'};
for i1 = 1:length(lib)
    if isempty(strfind(path, lib{i1}));
        path(lib{i1},path);
    end
end

% Directory name segment for Existing realisations/ensemble members -- just
% provide manually 
dns_r = {... cmip6
    '1001.001'; '1021.002'; '1041.003'; '1061.004'; '1081.005'; ...
    '1101.006'; '1121.007'; '1141.008'; '1161.009'; '1181.010'; ...
    '1231.001'; '1231.002'; '1231.003'; '1231.004'; '1231.005'; ...
    '1231.006'; '1231.007'; '1231.008'; '1231.009'; '1231.010'; ...
    '1251.001'; '1251.002'; '1251.003'; '1251.004'; '1251.005'; ...
    '1251.006'; '1251.007'; '1251.008'; '1251.009'; '1251.010'; ...
    '1281.001'; '1281.002'; '1281.003'; '1281.004'; '1281.005'; ...
    '1281.006'; '1281.007'; '1281.008'; '1281.009'; '1281.010'; ...
    '1301.001'; '1301.002'; '1301.003'; '1301.004'; '1301.005'; ...
    '1301.006'; '1301.007'; '1301.008'; '1301.009'; '1301.010'; ...
    ... 04.05.21. the rest of the 50 with the OTHER FORCING smbb
    '1011.001'; '1031.002'; '1051.003'; '1071.004'; '1091.005'; ...
    '1111.006'; '1131.007'; '1151.008'; '1171.009'; '1191.010'; ...
    '1231.011'; '1231.012'; '1231.013'; '1231.014'; '1231.015'; ...
    '1231.016'; '1231.017'; '1231.018'; '1231.019'; '1231.020'; ...
    '1251.011'; '1251.012'; '1251.013'; '1251.014'; '1251.015'; ...
    '1251.016'; '1251.017'; '1251.018'; '1251.019'; '1251.020'; ...
    '1281.011'; '1281.012'; '1281.013'; '1281.014'; '1281.015'; ...
    '1281.016'; '1281.017'; '1281.018'; '1281.019'; '1281.020'; ...
    '1301.011'; '1301.012'; '1301.013'; '1301.014'; '1301.015'; ...
    '1301.016'; '1301.017'; '1301.018'; '1301.019'; '1301.020'}; 
nr = length(dns_r);

% File name segment wrt. time - historical, same for each realisation
fns_th = {...
    '185001-185912'; '186001-186912'; '187001-187912'; '188001-188912'; ...
    '189001-189912'; '190001-190912'; '191001-191912'; '192001-192912'; ... 
    '193001-193912'; '194001-194912'; '195001-195912'; '196001-196912'; ...
    '197001-197912'; '198001-198912'; '199001-199912'; '200001-200912'; ...
    '201001-201412'};
ntch = length(fns_th); % number of time series chunks - historical

% same but for the SSP370 scenario
fns_ts = {...
    '201501-202412'; '202501-203412'; '203501-204412'; '204501-205412'; ...
    '205501-206412'; '206501-207412'; '207501-208412'; '208501-209412'; ...
    '209501-210012'};
ntcs = length(fns_ts); % * - scenario

ny = 2100-1850+1; 

% 05.04.19. First off, get the dimensions of the data, just opening one
% file.
% Deleted stuff for SST here due to redundancy

% Same thing about precipitation 03.05.21. which will do for surf
% temperature too 0.05.21. and also SST regridded onto ocean (variable SST
% not TEMP)
file1 = ['/proj/jedwards/archive/b.e21.BHISTcmip6.f09_g17.LE2-1001.001/atm/proc/tseries/month_1/' ...
    'b.e21.BHISTcmip6.f09_g17.LE2-1001.001.cam.h0.PRECT.185001-185912.nc'];
% For the atmosphere the grid is regular - Gaussian, i guess -- and so lat,
% lon are 1D var's. 
lat1 = ncread(file1,'lat');
lon1 = ncread(file1,'lon');   
% 19.01.20.
nlat1 = length(lat1);
nlon1 = length(lon1);

% 04.05.21. Define the Nino34 box in terms of lat1 instead of lat0
latmin = -5; latmax = 5; lonmin = 360-170; lonmax = 360-120; % Nino34
i_lat3 = find(lat1>latmin & lat1<latmax);
i_lon3 = find(lon1>lonmin & lon1<lonmax);
[lon3,lat3] = meshgrid(lon1(i_lon3),lat1(i_lat3)); 
    % lat is first dimension, unlike for ncread -> transpose e.g. at areal
    % averaging

% 05.04.19. New loop on realisations
nino34_djf = zeros(ny,nr);
nino34_jja = zeros(ny,nr);
PR_djf = zeros(nlon1,nlat1,ny,nr);
PR_jja = zeros(nlon1,nlat1,ny,nr);
TS_djf = zeros(nlon1,nlat1,ny,nr);
TS_jja = zeros(nlon1,nlat1,ny,nr);
parfor i1 = 1:nr
%for i1 = 1:nr
    'loop on realisation'
    i1
    
    % SST
    % 04.05.21. Everything i need is in the atm directories. But there are
    % two sets belonging to the two forcings. 
    if i1 < 51
    forcing = 'cmip6';
    else
    forcing = 'smbb';
    end
    datadir = ['/proj/jedwards/archive/b.e21.BHIST' forcing '.f09_g17.LE2-' dns_r{i1} '/atm/proc/tseries/month_1/'];
    
    % Loop for time chunks
    xh = [];
    for i2 = 1:ntch
        %i2
        file = [datadir 'b.e21.BHIST' forcing '.f09_g17.LE2-' dns_r{i1} '.cam.h0.SST.' fns_th{i2} '.nc'];
        xi2  = squeeze(ncread(file,'SST'));
        xh = cat(3,xh,xi2);
    end
    
    datadir = ['/proj/jedwards/archive/b.e21.BSSP370' forcing '.f09_g17.LE2-' dns_r{i1} '/atm/proc/tseries/month_1/'];
   
    xs = [];
    for i2 = 1:ntcs
        file = [datadir 'b.e21.BSSP370' forcing '.f09_g17.LE2-' dns_r{i1} '.cam.h0.SST.' fns_ts{i2} '.nc'];
        xi2  = squeeze(ncread(file,'SST'));
        xs = cat(3,xs,xi2); % concatenate wrt. time
    end
    x = cat(3,xh,xs); 
    
    % 04.05.21. 
    x_mm = movmean(x,3,3);

    nino34_djf(:,i1) = nansum(x_mm(i_lon3,i_lat3,1:12:end).*repmat(cosd(lat3'),1,1,ny),[1,2])/sum(cosd(lat3),[1,2]);
    nino34_jja(:,i1) = nansum(x_mm(i_lon3,i_lat3,7:12:end).*repmat(cosd(lat3'),1,1,ny),[1,2])/sum(cosd(lat3),[1,2]);
    
    % Precip
    datadir = ['/proj/jedwards/archive/b.e21.BHIST' forcing '.f09_g17.LE2-' dns_r{i1} '/atm/proc/tseries/month_1/'];
    
    xh = [];
    for i2 = 1:ntch
        file = [datadir 'b.e21.BHIST' forcing '.f09_g17.LE2-' dns_r{i1} '.cam.h0.PRECT.' fns_th{i2} '.nc'];
        xi2  = ncread(file,'PRECT');
        xh = cat(3,xh,xi2);
    end
    
    datadir = ['/proj/jedwards/archive/b.e21.BSSP370' forcing '.f09_g17.LE2-' dns_r{i1} '/atm/proc/tseries/month_1/'];
   
    xs = [];
    for i2 = 1:ntcs
        file = [datadir 'b.e21.BSSP370' forcing '.f09_g17.LE2-' dns_r{i1} '.cam.h0.PRECT.' fns_ts{i2} '.nc'];
        xi2  = ncread(file,'PRECT');
        xs = cat(3,xs,xi2); % concatenate wrt. time
    end
    y = cat(3,xh,xs); 
    
    % 03.05.21. 
    y_mm = movmean(y,3,3);
    PR_djf(:,:,:,i1) = y_mm(:,:,1:12:end);
    PR_jja(:,:,:,i1) = y_mm(:,:,7:12:end);
    
    % 03.05.21. Surface temperature
    datadir = ['/proj/jedwards/archive/b.e21.BHIST' forcing '.f09_g17.LE2-' dns_r{i1} '/atm/proc/tseries/month_1/'];
    
    xh = [];
    for i2 = 1:ntch
        file = [datadir 'b.e21.BHIST' forcing '.f09_g17.LE2-' dns_r{i1} '.cam.h0.TS.' fns_th{i2} '.nc'];
        xi2  = ncread(file,'TS');
        xh = cat(3,xh,xi2);
    end
    
    datadir = ['/proj/jedwards/archive/b.e21.BSSP370' forcing '.f09_g17.LE2-' dns_r{i1} '/atm/proc/tseries/month_1/'];
   
    xs = [];
    for i2 = 1:ntcs
        file = [datadir 'b.e21.BSSP370' forcing '.f09_g17.LE2-' dns_r{i1} '.cam.h0.TS.' fns_ts{i2} '.nc'];
        xi2  = ncread(file,'TS');
        xs = cat(3,xs,xi2); % concatenate wrt. time
    end
    y = cat(3,xh,xs); 
    
    % 03.05.21. 
    y_mm = movmean(y,3,3);
    TS_djf(:,:,:,i1) = y_mm(:,:,1:12:end);
    TS_jja(:,:,:,i1) = y_mm(:,:,7:12:end);
    
    % Uncomment this if 'for' used instead of 'parfor'
    %clear x xlin y %ylin % some large var's unnecessary to save
end
save('ENSO_IM_teleco_46_x.mat','-v7.3') 