#!/usr/bin/octave -qf

% this matlab/octave script will calculate the difference between two
% vibrationally resolved electronic spectra
%   Copyright 2015 by Jan von Cosel

warning('off', 'all');
pkg load signal;
format long;

arg_list = argv();
filename = arg_list{1}; % the input file MUST be the 1st argument
triangle_width = str2num(arg_list{2}); % the width of the triangular weighting function

% read in the spectrum
input_data = load(filename);

% extract the data
energies = input_data(:,1);
cold_spec = input_data(:,2);
hot_spec = input_data(:,3);

% create the x-axis for the correlation functions
span = length(energies) -1;   % the amount by which the correlation shifts the spectra (the "lag")
de = (energies(end) - energies(1)) / (span);
energies_long = (-((span) * de) : de : ((span) * de))';

% create the triangular weighting function
center = span + 1;    % the center point of the weighting function
triangle_weight = zeros(length(energies_long), 1);
for i = 1:length(energies_long)
    if abs(i - center) * de < triangle_width
        triangle_weight(i) = 1 - ((abs(i - center) * de) / triangle_width);
    endif
endfor

% create the simple difference weighting functions
simdiff_cross_weight = zeros(length(energies_long), 1);
simdiff_cross_weight(center) = 1;
simdiff_auto_weight = ones(length(energies_long), 1);

% calculate the correlation functions
auto_cold = xcorr(cold_spec) * de;
auto_hot = xcorr(hot_spec) * de;
cross = xcorr(cold_spec, hot_spec) * de;

% normalize the correlation functions
auto_cold_norm = auto_cold / (trapz(energies, cold_spec) * trapz(energies, cold_spec));
auto_hot_norm = auto_hot / (trapz(energies, hot_spec) * trapz(energies, hot_spec));
cross_norm = cross / (trapz(energies, cold_spec) * trapz(energies, hot_spec));

% calculate the similarity measure
S_cc = trapz(energies_long, triangle_weight .* auto_cold) ...
        / sqrt(trapz(energies_long, triangle_weight .* auto_cold) * trapz(energies_long, triangle_weight .* auto_cold));
S_hh = trapz(energies_long, triangle_weight .* auto_hot) ...
        / sqrt(trapz(energies_long, triangle_weight .* auto_hot) * trapz(energies_long, triangle_weight .* auto_hot));
S_ch = trapz(energies_long, triangle_weight .* cross) ...
        / sqrt(trapz(energies_long, triangle_weight .* auto_cold) * trapz(energies_long, triangle_weight .* auto_hot));

% calculate the difference
d_triangle = S_cc + S_hh - 2 * S_ch;
d_simdiff = auto_cold_norm(center) + auto_hot_norm(center) - 2 * cross_norm(center);

printf("%20.10f %20.10f\n", d_triangle, d_simdiff);
