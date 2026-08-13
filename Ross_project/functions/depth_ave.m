function [atten_depth,total_thickness]=depth_ave(atten_rate,depth,index)
% This function calculate the depth averaged attenuation at each depth
%           [atten_depth,total_thickness]=depth_ave(atten_rate,depth,index)
% Inputs:                                                                    
%           atten_rate = one-way attenuation at each depth (dB/km)
%           depth = time-converted depth (m) where surface = 0 m and negative downward
%           index = index of calculated section in boolean (1: within the range & 0: out of the range)
% Outputs:
%           atten_depth = depth-averaged attenuation rate at each depth (dB/km)
%               ** Note: values outside of ice thickness is set to NaN **
%           total_thickness = total ice thickness over the calculated section (m)

[rows,cols]=size(depth);
if rows<cols
    depth=depth';
end
atten_depth=NaN*ones(size(atten_rate));
thickness_section=depth(index); % all indices that we calculate
total_thickness=abs(thickness_section(end)-thickness_section(1));
dz=diff(depth);
output=cumsum(atten_rate(index).*abs(dz(index)),'omitnan')/total_thickness;
atten_depth(index)=output;

end