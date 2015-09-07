#!/usr/bin/octave -qf

% Copyright 2015 Jan von Cosel
%
% This file is part of utility-scripts.
%
% utility-scripts is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% utility-scripts is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have recieved a copy of the GNU General Public License
% along with utility-scripts. If not, see <http://www.gnu.org/licenses/>.


format long;

fs2au = 41.34137333656; % conversion between femtoseconds and atomic units of time
eh2cm = 4.55633526e-6;  % conversion between Hartree and cm-1

% read the data.
arg_list = argv();
filename = arg_list{1};
tau = str2num(arg_list{2}) * fs2au;
input_data = load(filename);

% find out the next largest power of two that we need for our spectrum.
N = 1;
while N < length(input_data)
  N = N * 2;
end

% construct the various times and frequencies needed.
dt = (input_data(2,1) - input_data(1,1)) * fs2au;
nyquist = 1.0 / (2.0 * dt);
time = (0.0 : dt : (N - 1) * dt)';
df = 2.0 * nyquist / N;
freq = (-nyquist : df : nyquist)';
wavenumbers = freq ./ (eh2cm / pi);

% construct the autocorrelation function with exponential damping.
corr = zeros(N,1);
corr(1:length(input_data)) = input_data(:,2) + j * input_data(:,3);
corr_damp = corr .* exp(-time / tau);

% this is the actual fourier transformation.
spec = fft(corr_damp);

% properly arrange the positive and negative parts of the spectrum.
spec_shift = zeros(N+1,1);
for i=0:(length(spec) / 2)
  spec_shift(i+1) = spec((length(spec) / 2)+i);
endfor
for i=1:(length(spec) / 2)
  spec_shift((length(spec) / 2)+1+i) = spec(i);
endfor

% write the results to the output files.
specfile = fopen("spec_out.dat", "w");
autofile = fopen("auto_out.dat", "w");
fprintf(autofile, "#         Time [au]        Re(Auto)         Im(Auto)\n");
fprintf(specfile, "# Frequency [au]    Frequency [cm-1]        Re(Spec)          Im(Spec)\n");
for i=1:N
  fprintf(autofile, "%20.10f %16.10f %16.10f\n", time(i), real(corr_damp(i)), imag(corr_damp(i)));
  fprintf(specfile, "%15.10f %20.10f %17.10f %17.10f\n", freq(i), wavenumbers(i), real(spec_shift(i)), imag(spec_shift(i)));
endfor
fprintf(specfile, "%15.10f %20.10f %17.10f %17.10f\n", freq(N+1), wavenumbers(N+1), real(spec_shift(N+1)), imag(spec_shift(N+1)));
fclose(specfile);
fclose(autofile);

