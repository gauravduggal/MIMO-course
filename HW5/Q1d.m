close all
clear all
clc

%sampling frequency
fs = 100e6;
%sampling time
ts = 1/fs;
%number of channel taps
N = 128;
%delay spread
tau_rms = 1000*1e-9;
%frequency sacing
delta_f = fs/N;
%delay axis
tau_axis = 0:ts:N*ts-ts;
%frequency_axis
f_axis = -fs/2:delta_f:fs/2-delta_f;

h = get_freq_selective_channel(tau_rms,fs,N);
%power
sum(abs(h).^2)


figure;
H = 1/(sqrt(N))*fftshift(fft(h,N));
legend_str ="Simulation: Delay spread:"+string(floor(tau_rms*1e9))+"ns";
semilogy(f_axis,abs(H),'-x','MarkerIndices',1:N/16:N,"DisplayName",legend_str,LineWidth=1.5)
xlabel("f (Hz) ")
ylabel("|H(f)|")
title("Frequency Domain channel amplitude ")
legend

figure;
plot(180/pi*angle(H),'-x','MarkerIndices',1:N/16:N,"DisplayName",legend_str,LineWidth=1.5)
xlabel("f (Hz)")
ylabel("phase (degree)")
title("Phase response")
title("Frequency domain channel phase")
legend

figure;
legend_str ="Simulation: Delay spread:"+string(floor(tau_rms*1e9))+"ns";
semilogy(tau_axis,abs(h).^2,'-^','MarkerIndices',1:N/16:N,"DisplayName",legend_str,LineWidth=1.5)
title("Power Delay Profile")
hold on

grid on

legend_str ="Theoretical: Delay spread:"+string(floor(tau_rms*1e9))+"ns";
gamma = 1/tau_rms;
semilogy(tau_axis,ts*gamma*exp(-tau_axis*gamma),'-o','MarkerIndices',1:N/16:N,"DisplayName",legend_str,LineWidth=1.5)
xlabel("\tau (ns)")
ylabel("|h(\tau)|^2")
legend


figure;
Ntrials = 10000;
val = zeros(1,Ntrials);
freq_corr = zeros(1,N);
for i=1:Ntrials
    h = get_freq_selective_channel(tau_rms,fs,N);
    H = 1/(sqrt(N))*fftshift(fft(h,N));
    val(i) = abs(H(N/2));
    freq_corr(1,:) = freq_corr(1,:) + 1/Ntrials*N*real(H(1)*conj(H));
%     sum(abs(h).^2)
end
legend_str ="Simulation: Delay spread:"+string(floor(tau_rms*1e9))+"ns";
plot(fs/2+f_axis,freq_corr,'-o','MarkerIndices',1:N/32:N,"DisplayName",legend_str,LineWidth=1.5)
xlabel("Frequency spacing (Hz)")
ylabel("Correlation")
title("Correlation vs Frequency spacing")
theory_axis = [linspace(0,fs/2,N/2),linspace(fs/2,0,N/2)]; 
freq_corr_theory = fft(ts*gamma*exp(-tau_axis*gamma),N);
% freq_corr_theory = abs(1./(1+1j*2*pi*theory_axis*tau_rms)); 
hold on
legend_str ="Theory: Delay spread:"+string(floor(tau_rms*1e9))+"ns";
plot(fs/2+f_axis,real(freq_corr_theory),'-x','MarkerIndices',1:N/32:N,"DisplayName",legend_str,LineWidth=1.5)
legend
grid on

figure;
histogram(val,Normalization="pdf")
histfit(val,30,'rayleigh')
legend("Simulation", "Theory")
title("Hist. of ampl. of a f bin of a freq. sel. chan.")
pd = fitdist(val','rayleigh')
xlabel("amplitude")
ylabel("frequency of occurence")
% histfit()

