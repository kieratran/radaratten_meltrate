function [geo_power, depth] = geo_correction(raw_power, time, surface_time, elevation)
% This function corrects returned power affected by geometric spreading
% Inputs:
%           raw_power = two-way returned power (Watts)
%           time = two-way time travel (s)
%           surface_time = two-way time return of the surface (s)
%           elevation = elevation of the aircraft (m)
% Output:
%           geo_power = two-way geometrically corrected power (dB)
%           depth = time-converted depth (m) where surface = 0 m and negative downward

v_ice = 1.68e8; % wave velocity through ice (m/s)
e_ice = 3.12; % ice permittivity

if isreal(raw_power) % convert from Watts to dB 
    converted_power = 10*log10(raw_power);
elseif imag(raw_power) ~= 0
    converted_power = 20*log10(abs(raw_power));
end

depth = -0.5*(time - surface_time)*v_ice;
geo_power = converted_power + 20*log10(2*(elevation+abs(depth)/sqrt(e_ice)));

end
