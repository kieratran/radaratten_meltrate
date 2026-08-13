function [attenuation] = temp2atten(temperature, E1, E2, E3, sig1, mu2, mu3, c2, c3)
    % Adapt from MacGregor 2007: Modeling Radar Attenuation at Siple Dome
    % Using Siple Dome constants for chemical impurities
    %
    % Inputs:
    %       Temperature in K
    %       %%%%%% Impurities based on MacGregor et al 2007
    %       num=1 is pure ice, num=2 is acidity, num=3 is sea salt
    %       E1; E2; E3; % Activation energy in eV
    %       sig1; % pure conductivity with 1M concentration (c1), unit uS/m
    %       mu2; mu3; % pure conductivity in S/m/M
    %       c2; c3; % radar-effective concentrations in uM (molarity)
    % Outputs:
    %       Attenuation through depth in dB/km
    %
    % [attenuation]=MacGregor(temperature)

    Tr = 251; % reference temp is -22.15 degree C
    k_bolt = 8.617333262e-5; % Boltzmann constant in eV/K
    esp0 = 8.85418782e-6; % permittivity of free space
    esp_prime = 3.12; % real part of permittivity of ice
    c = 299792458; % speed of light in m/s
    e = 2.718281828;

    % Do work
    temp_ratio = (1/Tr - 1./temperature)./k_bolt;

    pure = sig1 .* exp(E1 .* temp_ratio);
    acid = mu2 * c2 .* exp(E2 .* temp_ratio);
    salt = mu3 * c3 .* exp(E3 .* temp_ratio);

    sigma = pure + acid + salt;
    const = 1000*10*log10(e)/(esp0*c*sqrt(esp_prime));
    attenuation = const .* sigma;
end