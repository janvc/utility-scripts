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
c0 = 137.0359998;       % the speed of light in atomic units

% read the data.
arg_list = argv();
filename = arg_list{1};
tau = str2num(arg_list{2}) * fs2au;
offset = str2num(arg_list{3});
input_data = load(filename);

% find out the next largest power of two that we need for our spectrum.
N = 1;
while N < length(input_data)
  N = N * 2;
end

% construct the various times and frequencies needed.
dt = (input_data(2,1) - input_data(1,1)) * fs2au;
nyquist = 2.0 * pi / (2.0 * dt);
time = (0.0 : dt : (N - 1) * dt)';
double_time = (-(N-1)*dt : dt : N * dt)';
df = 2.0 * nyquist / N;
double_df = 2.0 * nyquist / (2*N);
freq = (-nyquist : df : nyquist)';
double_freq = (-nyquist : double_df : nyquist)';

% construct the autocorrelation function with cosine weighting and exponential damping.
corr = zeros(N,1);
corr(1:length(input_data)) = input_data(:,2) + j * input_data(:,3);
corr_damp = corr .* exp(-time / tau) .* cos(pi * time / (2.0 * time(N)));

% construct the mirrored autocorrelation function with cosine weighting and exponential damping
double_corr = zeros(2*N,1);
for a=0:N-1
  double_corr(a+1) = conj(corr(N-a));   % complex conjugate the negative part!!!
  double_corr(2*N-1-a) = corr(N-a);
end
double_corr_damp = double_corr .* cos(pi * double_time / (2.0 * time(N))) .* exp(-abs(double_time) / tau);

% this is the actual fourier transformation.
spec = ifft(corr_damp);
double_spec = ifft(double_corr_damp);

% properly arrange the positive and negative parts of the spectrum.
spec_shift = zeros(N+1,1);
double_spec_shift = zeros(2*N+1,1);
spec_shift(1:N/2) = spec((N/2)+1:N);
spec_shift((N/2)+1:N+1) = spec(1:(N/2)+1);
double_spec_shift(1:N) = double_spec(N+1:2*N);
double_spec_shift(N+1:2*N+1) = double_spec(1:N+1);

% write the results to the output files.
specfile = fopen("spec_out.dat", "w");
autofile = fopen("auto_out.dat", "w");
double_specfile = fopen("spec_out_double.dat", "w");
double_autofile = fopen("auto_out_double.dat", "w");
offset_specfile = fopen("spec_out_offset.dat", "w");
fprintf(autofile, "#         Time [au]        Re(Auto)         Im(Auto)\n");
fprintf(specfile, "# Frequency [au]    Frequency [cm-1]        Re(Spec)          Im(Spec)\n");
fprintf(double_autofile, "#         Time [au]        Re(Auto)         Im(Auto)\n");
fprintf(double_specfile, "# Frequency [au]    Frequency [cm-1]        Re(Spec)          Im(Spec)\n");
fprintf(offset_specfile, "# Frequency [au]        Re(Spec)          Im(Spec)\n");
for a=1:length(time)
  fprintf(autofile, "%20.10f %16.10f %16.10f\n", time(a), real(corr_damp(a)), imag(corr_damp(a)));
end
for a=1:length(freq)
  fprintf(specfile, "%15.10f %17.10f %17.10f\n", freq(a), real(spec_shift(a)), imag(spec_shift(a)));
end
for a=1:length(double_time)
  fprintf(double_autofile, "%20.10f %16.10f %16.10f\n", double_time(a), real(double_corr_damp(a)), imag(double_corr_damp(a)));
end
factor = 4.0 * pi^2 / (3.0 * c0);
for a=1:length(double_freq)
  fprintf(double_specfile, "%15.10f %17.10f %17.10f\n", double_freq(a), real(double_spec_shift(a)), imag(double_spec_shift(a)));
  fprintf(offset_specfile, "%15.10f %17.10f %17.10f %17.10f\n", double_freq(a)+offset, real(double_spec_shift(a))*factor*(double_freq(a)+offset), imag(double_spec_shift(a))*factor*(double_freq(a)+offset), abs(double_spec_shift(a))*factor*(double_freq(a)+offset));
end
fclose(specfile);
fclose(autofile);
fclose(double_specfile);
fclose(offset_specfile);
fclose(double_autofile);

