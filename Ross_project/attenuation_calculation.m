% ==== Author: Kiera Tran ====
% This script is used for calculating ice-penetrating radar attenuation
clear all; clc;

addpath(genpath('~/Ross_projects'))   % add general directory to all existing functions in Ross_project folder

%% Load data:
% In this section, you can load radargram add any directory or folder and load data
% Best way to test run the code is to use Operation IceBridge data: https://data.cresis.ku.edu/data/rds/
%     |---- Note: use CSARP_standard

% For example:
dataFile = 'radar_data.mat'; % Specify your data file here
load(dataFile, 'Data', 'surfLayer', 'bedLayer', 'time', 'lat', 'lon', 'plane_elevation', 'surfTime');   % Load necessary variables
% Surface and Bed layers can either be in elevation, index, or time
% In this case, 'surfLayer' and 'bedLayer' are indices of surface and base layers
%     |---- Check units and conversions
%         |---- elevation = -0.5 * (time - surface_time) * 1.68e8;   %  Sound velocity through ice is 1.68e8 m/s


%% Find inflection points:
% This section identify (1) the boundary between firn and ice, 
% and (2) bed-echo length to neglect scattering power signals

% ---- 1. Initialization ----
inLayer = nan(size(lat));   % indices of firn depth
inLayer2 = nan(size(lat));   % indices of bed-echo return
firn_confidence = nan(size(lat));

% ---- 2. Setting limits for each layers (you can tune this) ----
firn_prior = 70;   echo_prior = 50;   % expected depths (m)
firn_range = 30;   echo_range = 50;   % search ± range (m)
disp('>>>> Getting layers ...')
for i = 1:numel(lat)
    if isnan(surfLayer(i)) || isnan(bedLayer(i)) || bedLayer(i) <= surfLayer(i) + 10
        continue   % all outputs already NaN from pre-allocation
    end
    [power, depth] = geo_correction(Data(:,i), time(:,i), surfTime(i), plane_elevation(i)); 

% ---- 3.  Dectecting layers' boundaries ----
    surf_s = surfLayer(i);
    bed_s = bedLayer(i);
    depth_col = depth(surf_s:bed_s);
    power_col = power(surf_s:bed_s);
    n_samp = numel(depth_col);
    dz = abs(mean(diff(depth_col)));
    % --- Firn/ice: search from surface downward ---
    firn_lo = max(2, round((firn_prior - firn_range) / dz));
    firn_hi = min(n_samp-1, round((firn_prior + firn_range) / dz));
    if firn_hi > firn_lo + 2
        [bp, firn_confidence(i)] = piecewise_fit(depth_col, power_col, firn_lo, firn_hi, "firn");
        inLayer(i) = surf_s + bp - 1;
    else
        warning('nn=%d: firn window too narrow (dz=%.2fm)', i, dz);
    end
    % --- Ice/saline: search from bed upward ---
    echo_lo = max(2, n_samp - round((echo_prior + echo_range) / dz));
    echo_hi = min(n_samp-1, n_samp - round((echo_prior - echo_range) / dz));
    if echo_hi > echo_lo + 2
        % Flip segment so "surface side" is always segment 1 in the fit
        d_echo = depth_col(echo_lo:echo_hi);
        p_echo = power_col(echo_lo:echo_hi);
        [bp, ~] = piecewise_fit(d_echo, p_echo, 2, numel(d_echo) - 1, "saline");
        % Map back — bp is relative to saline segment
        inLayer2(i) = surf_s + echo_lo + bp - 2;
    else
        warning('nn=%d: saline window too narrow (dz=%.2fm)', i, dz);
    end
end

% ---- 4. Smoothing layers ----
nanidx = isnan(inLayer); inLayer = smoothn(inLayer, 600); 
inLayer(nanidx) = NaN; inLayer = round(inLayer);
nanidx = isnan(inLayer2); inLayer2 = smoothn(inLayer2, 600); 
inLayer2(nanidx) = NaN; inLayer2 = round(inLayer2);

% ---- 5. Converting from indices to elevation (m) ----
firnElevation = nan(size(lat));
echoElevation = nan(size(lat));
surfElevation = nan(size(lat));
bedElevation = nan(size(lat));
for i = 1:numel(lat)
    if ~isnan(surfTime(i)) && ~isnan(inLayer(i)) && ~isnan(inLayer2(i))
        [power, depth] = geo_correction(Data(:,i), time(:,i), surfTime(i), plane_elevation(i));
        firnElevation(i) = depth(inLayer(i));   % firn layer depth (m)
        echoElevation(i) = depth(inLayer2(i));   % bed-echo layer depth (m)
        surfElevation(i) = depth(surfLayer(i));   % surface layer depth (m)
        bedElevation(i) = depth(bedLayer(i));    % bed layer depth (m)
    end
end

% ---- 6. Calculating ice thickness (m) ----
IceThick = abs(bedElevation - surfElevation);

%% Calculate attenuation rate:
% This section calculates depth-resolved multi-reflector attenuation rates

% ---- 1. Initialization ----
AttenRate = nan(size(lat));   % linear fitting attenuation rate
AttenRate_unc = nan(size(lat));   % linear fitting uncertainty
Na = nan(size(lat));   % piecewise fitting attenuation rate
Na_unc = nan(size(lat));    % piecewise fitting uncertainty

disp('>>> Calculating attenuation rates ...')
for i = 1:numel(lat)
    if ~isnan(inLayer(i)) && ~isnan(inLayer2(i)) && inLayer2(i) > inLayer(i)   % surface and bed picks exist
        [power, depth] = geo_correction(Data(:,i), time(:,i), surfTime(i), plane_elevation(i));   % geometrically correct returned power

% ---- 2. Looping through each segment length and calculate attenuation ----
        num = 1;
        for ii = 20:5:50 % segment length is based on percentage of total ice thickness (20% to 50% of ice thickness)
            % --- Piecewise attenuation rate
            atten_rate = atten_calc(power, depth, inLayer(i), inLayer2(i), 0, ii);
            atten_depth = depth_ave(atten_rate, depth, inLayer(i):inLayer2(i));
            Na(num, i) = atten_depth(inLayer2(i));
            % --- Piecewise uncertainty
            atten_rate = piecewise_unc(power, depth, inLayer(i), inLayer2(i), 0, ii);
            atten_depth = depth_ave(atten_rate, depth, inLayer(i):inLayer2(i));
            Na_unc(num, i) = atten_depth(inLayer2(i));
            num = num+1;
        end
        % --- Segment length = 100% ice thickness & uncertainty
        p = polyfit(depth(inLayer(i):inLayer2(i)), power(inLayer(i):inLayer2(i)),1)/2*1000;   % One-way linear fitting attenuation rate (dB/km)
        AttenRate(i) = p(1,1);
        AttenRate_unc(i) = slopeSE(depth(inLayer(i):inLayer2(i)), power(inLayer(i):inLayer2(i)));   % uncertainty
        % --- Combine all segment lengths
        Na(num, i) = AttenRate(i);
        Na_unc(num, i) = AttenRate_unc(i);
    end
end

% ---- 3. Averaging depth-resolved attenuation rates with different piecewise regressions
 Attenuation = nanmean(Na);
 Attenuation_unc = nanmean(Na_unc);


%% Sanity check and visualization:
% This will load A-scope (radargram), 
% and you can choose specific Z-scope (power profile) to view for validation
plot_radargram