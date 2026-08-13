function atten_rate = atten_calc(geo_power, depth, surface_pick, bed_pick, segment_length, percentage)
% This function calculate the attenuation at each depth within the profile
% Inputs:
%           geo_power = two-way geometrically corrected power (dB)
%           depth = time-converted depth (m) where 0 is at surface and negative is downword
%           surface_pick = index of surface pick respective to time travel
%           bed_pick = index of bottom pick respective to time travel
%           segment_length = length of each segment (m)
%           percentage = percent of ice thickness (%)
% Output:
%           atten_rate = one-way attenuation at each depth (dB/km)
%               ** Note: values outside of ice thickness is set to NaN **

atten_rate = NaN*ones(size(geo_power));

if segment_length == 0 && percentage ~= 0 % choose segment length based on percent of ice thickness
    segment_length = abs(depth(bed_pick) - depth(surface_pick)) * percentage/100;
elseif segment_length ~=0 && percentage ~=0
    error('Conflict arrived!! Make up your mind.')
elseif segment_length == 0 && percentage == 0
    error('Nothing can be run if there is no segment length')
end

if abs(depth(bed_pick) - depth(surface_pick))/segment_length >= 1 % ice thickness must be greater than segment length
    starting = depth(surface_pick);
    ending = starting - segment_length;
    n = 1;
    while ending >= depth(bed_pick) % keep sliding the window until segment is smaller than the window
        index = find(depth <= (starting + 0.00005) & depth >= (ending - 0.00005)); % +/- 0.00005 for expanding the range
        p = polyfit(depth(index), geo_power(index), 1);
        atten_rate(index) = p(1, 1);
        n = n+1;
        starting = ending;
        ending = starting - segment_length;
        %plot(polyval(p, depth(index)), depth(index), 'k', 'LineWidth', 1.5); hold on;
    end
    ending = depth(bed_pick); % running the leftover segment of the profile
    index = find(depth <= (starting + 0.00005) & depth >= (ending - 0.00005));
    p = polyfit(depth(index), geo_power(index), 1);
    atten_rate(index)=p(1,1);
    %plot(polyval(p, depth(index)), depth(index), 'k', 'LineWidth', 1.5, 'DisplayName', 'Piecewise'); hold on;
end
atten_rate = atten_rate / 2 * 1000; % convert from two- way dB/m to one-way dB/km



end