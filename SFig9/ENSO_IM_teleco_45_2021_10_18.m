% 03.05.21. Developed from *_23.m. Companion to *_22.m. Calculate the
% E-wise correlations etc. for each year. 

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

% 04.05.21. Load the variables nino34_djf nino34_jja PR_djf PR_jja TS_djf
% TS_jja one at a time ...
% switches
sw_s = 1; % 1: djf; 2: jja
sw_v = 2; % 1:  pr; 2:  ts

season = {'djf'; 'jja'};
prorts = {'PR'; 'TS'};
varname = [prorts{sw_v} '_' season{sw_s}];
ninoname = ['nino34_' season{sw_s}];
load('ENSO_IM_teleco_46.mat',varname,ninoname,'nlon1','nlat1','ny') % takes ~2 minutes
% ... then use a new variable of a generic name
eval(['myvar = '  varname ';'])
eval(['nino  = ' ninoname ';'])

% Data variables to save
% 25.03.21. Calculate the E-std so that i can evaluate my formula for the
% variance change explained by the change of variance explained by ENSO.
% This calculation can be outside of the loop as Matlab's std can act along
% a dim'n
sd = squeeze(std(myvar,[],4));
%
cc = zeros(nlon1,nlat1,ny); % corr coeff
rc = cc; % reg coeff
sr = cc; % std of residuals
for i1 = 1:nlon1
    i1
    for i2 = 1:nlat1
        % Make an inner loop a parfor because we don't want to broadcast
        % huge arrays like P
        cc_i2 = zeros(ny,1); % corr coeff
        rc_i2 = cc_i2; % reg coff
        sr_i2 = cc_i2; % std of residuals
        myvar_i2 = squeeze(myvar(i1,i2,:,:));
        parfor i3 = 1:ny
        %for i3 = 1:ny
            %i3
            x = nino(i3,:)';
            y = myvar_i2(i3,:)';
            temp = corrcoef(x,y);
            cc_i2(i3) = temp(1,2);
            % begin: Can skip these to save time or if not needed
            tbl = table(y,x);
            mdl = fitlm(tbl,'y ~ x');
            rc_i2(i3) = mdl.Coefficients{2,1};   
            sr_i2(i3) = mdl.MSE;
            % end: Can skip these ...
        end   
        cc(i1,i2,:,:) = cc_i2;
        rc(i1,i2,:,:) = rc_i2;
        sr(i1,i2,:,:) = sr_i2;
    end
end


save(['ENSO_IM_teleco_45_' prorts{sw_v} '_' season{sw_s} '_x.mat'],...
    '-v7.3','sd','cc','rc','sr','tw','rm')