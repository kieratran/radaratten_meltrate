% This script is also run in attenuation_calculation.m
% Use for plotting radargram 
disp('Converting radar power to dB ........')


disp('>>>> Plotting radargram ...')



f = figure('color','White');
imagesc(20*log10(abs(Data)));
colormap bone
colorbar
hold on;
plot(surfLayer,'b','LineWidth',1);
plot(bedLayer,'r','LineWidth',1);
plot(inLayer, 'y', 'LineWidth', 1);
plot(inLayer2, 'g', 'LineWidth', 1);
saveas(f,flightname,'jpg');
close all;


f = figure('color','White');
imagesc(20*log10(abs(Data))); hold on;   % if Data is complex
% imagesc(10*log10(abs(Data))); hold on;   % if Data is real 
colormap bone; colorbar; 
plot(surfLayer, 'b', 'LineWidth', 1); 
plot(bedLayer, 'r', 'LineWidth', 1);
saveas(f, flightname, 'jpg'); % save radargram image as JPEG

[power, depth] = geo_correction(Data(:, i), time(:, i), surfTime(i), elev(i));
% If geo_correction returns linear amplitude, uncomment:
% power = 20*log10(abs(power));

figure('Name', sprintf('A-scope — %s  nn=%d', flightname, nn), ...
       'NumberTitle', 'off', 'Position', [100 100 500 700]);

plot(power, depth, 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
hold on;

% Crop to surface–bed + 10% padding
d_surf = surfElevation(nn);
d_bed  = bedElevation(nn);
if ~isnan(d_surf) && ~isnan(d_bed)
    pad = abs(d_surf - d_bed) * 0.1;
    ylim([d_bed - pad, d_surf + pad]);
end

% --- Layer markers ---
layers = { ...
    surfElevation(nn),  'Surface',           [0.00 0.45 0.70]; ...
    bedElevation(nn),   'Bed',               [0.80 0.20 0.20]; ...
    inElevation(nn),    'Firn/ice',    [0.93 0.69 0.13]; ...
    inElevation2(nn),   'Ice/saline',  [0.47 0.67 0.19]; ...
};

for k = 1:size(layers, 1)
    e   = layers{k, 1};
    lbl = layers{k, 2};
    col = layers{k, 3};
    if ~isnan(e)
        yline(e, '--', 'Color', col, 'LineWidth', 1.6, ...
              'Label', lbl, ...
              'LabelVerticalAlignment', 'bottom', ...
              'LabelHorizontalAlignment', 'right', ...
              'FontSize', 9);
    end
end

xlabel('Power (dB)');
ylabel('Elevation (m)');
title(sprintf('%s   nn=%d   Lat %.4f  Lon %.4f', ...
      flightname, nn, Lat{n}(nn), Lon{n}(nn)), ...
      'Interpreter', 'none', 'FontSize', 10);
grid on;