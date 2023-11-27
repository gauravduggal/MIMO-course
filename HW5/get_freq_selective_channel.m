function [h] = get_freq_selective_channel(tau_rms,fs,N)
ts = 1/fs;
h = randn(1,N).*sqrt((0.5*(1-exp(-ts/tau_rms))/(1-exp(-N*ts/tau_rms)) * exp(-1*(0:N-1)*ts/tau_rms))) ... 
+ 1j*randn(1,N).*sqrt((0.5*(1-exp(-ts/tau_rms))/(1-exp(-N*ts/tau_rms)) * exp(-1*(0:N-1)*ts/tau_rms)));
end