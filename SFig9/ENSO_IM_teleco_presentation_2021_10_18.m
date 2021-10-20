% 01.02.21. Plot diagrams for the presentation paper. 

% 05.02.21. Write fields of corr coeff and its 21st c. slope into nc files
% so that one could make plots in panoply as an alternative to Matlab. 

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
lib =    {'/proj/bodai/MATLAB/Useful'; ...
          '/proj/bodai/MATLAB/Useful/chadagreene-CDT-39b51d9/cdt'};
for i1 = 1:length(lib)
    if isempty(strfind(path, lib{i1}))
        path(lib{i1},path);
    end
end

mkvid = 1;
% A path where movie frames as jpg files can be saved to and loaded from
if mkvid
path_4_frames = ['/proj/bodai/MATLAB/ENSO_IM_teleco_presentation_' datestr(now,30) '/'];
system(['mkdir ' path_4_frames]);
end

S_land = shaperead('landareas','UseGeoCoords',true);


load ENSO_IM_teleco_41_12_PR_jja
load ENSO_IM_teleco_46 lat1 lon1

fh = figure;
colormap(redblue)
% set(fh,'Position',[100 100 900 300]) 
% ax_pos = [... avoid overlaps, otherwise Matlab might not plot!!!
%     0.01 0.01 0.45 0.96;
%     0.5  0.01 0.45 0.96];
% 12.02.21. For a vertical stack
set(fh,'Position',[100 100 400 600]) 
ax_pos = [... avoid overlaps, otherwise Matlab might not plot!!!
    0.01 0.50 0.96 0.45;
    0.01 0.01 0.96 0.45];

rg_max = 1;

mytitle = {'mean sgn(r)r^2 in 1980-2014'; 'difference of r^2 bw. 1980-2014 and 2066-2100'};
    % 25.03.21. changed r -> r^2.
for i1 = 1:2
    figure(fh)

    subplot('Position',ax_pos(i1,:))
    ax_pos(i1,:)
    ax = axesm('Robinson','Origin',[0 180 0]);
    framem; gridm;
    axis off
    axis(ax,'tight')

    hcb = colorbar('southoutside');
    caxis(rg_max*[-1 1]);
    title(mytitle{i1});

    % Main plot command
    if i1 == 1
    %pcolorm(lat1, lon1, cc_ave_1'); 
    %pcolorm(lat1, lon1, cc_ave_1'.^2);
    pcolorm(lat1, lon1, cc_ave_1'.^2.*sign(cc_ave_1')); 
    else
    %pcolorm(lat1, lon1, cc_ave_2' - cc_ave_1');
    pcolorm(lat1, lon1, cc_ave_2'.^2 - cc_ave_1'.^2);
    %pcolorm(lat1, lon1, cc_ave_2'.^2.*sign(cc_ave_2') - cc_ave_1'.^2.*sign(cc_ave_1')); 
    end
        % 24.06.20. Origin set above; and no need to wrap around running
        % lon from 0 to 360

    % stiple for NON-significance!
    if i1 == 2
    stipplem(lat1,lon1,pt > 0.05,'MarkerSize',2,'Color','m'); 
    end

    geoshow([S_land.Lat ], [S_land.Lon],'Color','black');
end

if mkvid
% Print into files -- could be frames of a movie
set(gcf,'paperpositionmode','auto')
print(gcf,'-r300','-djpeg',[path_4_frames 'ENSO_IM_teleco_41']);
saveas(gcf,[path_4_frames 'ENSO_IM_teleco_41'])
end

% 05.02.21. Variables to be written out to nc file
nlon1 = length(lon1);
nlat1 = length(lat1);

ncfn = [path_4_frames 'ENSO_IM_teleco_presentation_2intervals.nc'];

nccreate(ncfn,'cc_ave_1','Dimensions',{'lon',nlon1,'lat',nlat1});
nccreate(ncfn,'cc_ave_2','Dimensions',{'lon',nlon1,'lat',nlat1});
nccreate(ncfn,'pt'      ,'Dimensions',{'lon',nlon1,'lat',nlat1});
nccreate(ncfn,'lon'     ,'Dimensions',{'lon',nlon1});                   
nccreate(ncfn,'lat'     ,'Dimensions',{'lat',nlat1});
nccreate(ncfn,'sd_ave_1','Dimensions',{'lon',nlon1,'lat',nlat1});
nccreate(ncfn,'sd_ave_2','Dimensions',{'lon',nlon1,'lat',nlat1});

ncwrite(ncfn,'cc_ave_1',cc_ave_1); 
ncwrite(ncfn,'cc_ave_2',cc_ave_2);
ncwrite(ncfn,'pt'      , pt);
ncwrite(ncfn,'lon'     ,lon1 );
ncwrite(ncfn,'lat'     ,lat1 );
ncwrite(ncfn,'sd_ave_1',sd_ave_1); 
ncwrite(ncfn,'sd_ave_2',sd_ave_2);

% 25.03.21. This is how the ENSO teleconnections can be related to variance
% changes
fh = figure;
colormap(redblue)
set(fh,'Position',[100 100 800 600]) 
ax_pos = [... avoid overlaps, otherwise Matlab might not plot!!!
    0.01 0.50 0.45 0.45;
    0.01 0.01 0.45 0.45;
    0.50 0.50 0.45 0.45;
    0.50 0.01 0.45 0.45];

var_1 = sd_ave_1'.^2; 
var_2 = sd_ave_2'.^2 - sd_ave_1'.^2;
%var_3 = ((sd_ave_2.^2.*cc_ave_2.^2 - sd_ave_1.^2.*cc_ave_1.^2) ./ var_2')';
 var_3 = ((sd_ave_2.^2.*cc_ave_2.^2 - sd_ave_1.^2.*cc_ave_1.^2)           )';
var_4 = (var_2' - (sd_ave_2.^2.*cc_ave_2.^2 - sd_ave_1.^2.*cc_ave_1.^2))';
%var_4 = var_3./var_4;
%rg_max_1 = max(abs(var_1(:))); % a location in indonesia stands out
rg_max_1 = quantile(    var_1(:) ,0.995);
rg_max_4 = quantile(abs(var_4(:)),0.995);

mytitle = {'mean var \sigma_1^2 in 1980-2014'; ...
    'difference \Delta\sigma^2=\sigma_2^2-\sigma_1^2 bw. 1980-2014 and 2066-2100'; ...
    '\Delta\sigma^2r^2=\sigma^2_2r^2_2-\sigma^2_1r^2_1'; '\Delta\sigma^2-\Delta\sigma^2r^2'};
for i1 = 1:4
    figure(fh)

    subplot('Position',ax_pos(i1,:))
    ax_pos(i1,:)
    ax = axesm('Robinson','Origin',[0 180 0]);
    framem; gridm;
    axis off
    axis(ax,'tight')

    hcb = colorbar('southoutside');
    title(mytitle{i1});

    % Main plot command
    if     i1 == 1
    caxis(rg_max_1*[-1 1]);
    pcolorm(lat1, lon1, var_1); 
    elseif i1 == 2
    caxis(rg_max_1*[-1 1]);
    pcolorm(lat1, lon1, var_2); 
    elseif i1 == 3 % see email 17.03.21.: fraction of variance change explained by change of variance explained by ENSO
    %caxis([-1 1]); actually, it can be well larger than 1 in modulus
    caxis(rg_max_1*[-1 1]);
    pcolorm(lat1, lon1, var_3);     
    else % in absolute terms: variance change unexpalined by *
    caxis(rg_max_1*[-1 1]);
    %caxis(rg_max_4*[-1 1]);
    pcolorm(lat1, lon1, var_4);    
    end
        % 24.06.20. Origin set above; and no need to wrap around running
        % lon from 0 to 360

    geoshow([S_land.Lat ], [S_land.Lon],'Color','black');
end

if mkvid
% Print into files -- could be frames of a movie
set(gcf,'paperpositionmode','auto')
print(gcf,'-r300','-djpeg',[path_4_frames 'ENSO_IM_teleco_41_teleco_n_std']);
saveas(gcf,[path_4_frames 'ENSO_IM_teleco_41_teleco_n_std'])
end