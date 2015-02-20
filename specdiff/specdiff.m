% this matlab/octave script will calculate the difference between two
% vibrationally resolved electronic spectra
%   Copyright 2015 by Jan von Cosel

format long;

% read in the raw spectra
raw_spectra = load("spectra-16.dat");
energies    = raw_spectra(:,1);
gsspec_raw  = raw_spectra(:,2);
esspec_raw  = raw_spectra(:,3);

% normalize the spectra
gsspec_norm = gsspec_raw / trapz(energies, gsspec_raw);
esspec_norm = esspec_raw / trapz(energies, esspec_raw);

% calculate the simple difference function
simple_diff = esspec_norm - gsspec_norm;


