function output=depth_ave_thickness(attenuation,thick)

output = NaN*ones(size(attenuation));
z = linspace(0,thick,length(attenuation));
dz = diff(z);
output(2:end) = cumsum(attenuation(2:end) .* abs(dz), 'omitnan')/thick;
output(1) = 0;

end
 