% 09.02.21. Companion to *_45, 46.m. Calculate the mean of the E-wise corr
% coeff in two periods and evaluate the significance whether the true neans
% are different performing a two-sample t-test. Need to use the Fisher
% transformation in order to satisfy the assumptions of the t-test.
% Reference for the methodology is:
% 
% Bódai, T., G. Drótos, M. Herein, F. Lunkeit, and V. Lucarini (2020) The
% Forced Response of the El Niño–Southern Oscillation–Indian Monsoon
% Teleconnection in Ensembles of Earth System Models. J. Climate, 33,
% 2163–2182, https://doi.org/10.1175/JCLI-D-19-0341.1
%
% There we performed a Mann-Kendall test, but the idea is the same: we want
% to test for nonstationarity, and for that the time seris is better made
% up of independent Gaussian random variables, differing at most wrt. their
% mean, i.e., no heteroscedasrticity is present.

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

% 04.05.21. 
varname = 'TS_jja';
load ENSO_IM_teleco_46 lat1 lon1 nlon1 nlat1 
load(['ENSO_IM_teleco_45_' varname])

% File name to create here or load for plotting. 09.02.21. Use some
% data set number manually input after 'ENSO_IM_teleco_41_' (12 here).
fname = ['ENSO_IM_teleco_41_12_' varname '.mat']; 

% Select the time periods over which we determine means and take
% differences
yr0 = 1850;
y10 = -yr0+1961;%1980; 
y1f = -yr0+1990;%2014; 
y20 = -yr0+2071;%2066; 
y2f = -yr0+2100; 

i2 = 1; % just use the zero lag data (18.10.21. In fact, that's all we have now from *_45.m)

cc_ave_1 = mean(cc(:,:,y10:y1f,i2),3);
cc_ave_2 = mean(cc(:,:,y20:y2f,i2),3);

% 25.03.21.
sd_ave_1 = mean(sd(:,:,y10:y1f),3);
sd_ave_2 = mean(sd(:,:,y20:y2f),3);

pt       = zeros(nlon1,nlat1);
for i4 = 1:nlon1 % took inside to make it parfor-loop
    %i2
    i4
    % Define these so that not the whole of cc is broadcast to parfor
    % workers. Might as well apply the Fisher transform here, needed
    % for the t-test.
    cc_1_i4 = atanh(squeeze(cc(i4,:,y10:y1f,i2))); 
    cc_2_i4 = atanh(squeeze(cc(i4,:,y20:y2f,i2)));
    pt_i4 = zeros(nlat1,1);
    %parfor i3 = 1:nlat1
    for i3 = 1:nlat1
        [~,p] = ttest2(cc_1_i4(i3,:),cc_2_i4(i3,:));
        pt_i4(i3) = p;
    end
    pt(i4,:) = pt_i4;
end

save(fname,'pt','y10','y1f','y20','y2f','cc_ave_1','cc_ave_2') 

save(fname,'sd_ave_1','sd_ave_2','-append') % 25.03.21.