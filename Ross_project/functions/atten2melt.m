function [tmd_Na_SD, tmd_Na_TD] = atten2melt(thickness, velocity, temp, smb, Ts, basal_slope, surface_slope, m)
% Parameters you do not need to change
dsdx = 0; % surface slope (only useful for non-ice-shelf; turned off)
dz = 5; % meter, vertical grid spacing
G = 0; % no geothermal heat flux needed for ice shelf; turned off

% Calibrations:
for ii = 1:length(m)
    tmd = thermal1d(thickness, dsdx, dz, smb + surface_slope, m(ii) + basal_slope, Ts, G, velocity, temp, "Nye", "FD");
    tmd = tmd.solveFI();
    solsol = tmd.sol.y(1,:);
    % Siple Dome:
    tmd_attenrate = temp2atten(solsol, 0.51, 0.2, 0.19, 9.2, 3.2, 0.43, 1.3, 4.2);
    tmd_attendepth = depth_ave_thickness(tmd_attenrate, thickness);
    tmd_Na_SD(ii) = tmd_attendepth(end);
    % Taylor Dome:
    tmd_attenrate = temp2atten(solsol, 0.51, 0.2, 0.19, 9.2, 3.2, 0.43, 1.57, 0.73); 
    tmd_attendepth = depth_ave_thickness(tmd_attenrate, thickness);
    tmd_Na_TD(ii) = tmd_attendepth(end);
end
%plot(b,tmd_Na,'k','LineWidth',2); hold on;

end