close all
clear all
clc

Pavg= 1;
fs = 1e6;


%fast time
ts = 1/fs;
T = 0.1;
t = 0:ts:T-ts;
%delay spread
tau_rms = 1e-6;
%Doppler spread
fdmax = 100;
%sub-carriers
N = 128;

%slow time

TS = N*ts;
slow_time = 0:TS:T-TS;

%Resolveable paths i.e. channel taps
L = 8;
Mr = 2;
Nt = 2;

lambda = 4;
d = 2;

%nonresolvable paths
N = 30;



Pt = 1;




phi = 360*rand(N,L);

AoA= 360*rand(1,N);
AoD = sqrt(360)*randn(1,N);
a_rx = exp(-1j*transpose(0:Mr-1)*2*pi*d/lambda*sind(AoA));
a_tx = transpose(exp(-1j*transpose((0:Nt-1))*2*pi*d/lambda*sind(AoD)));
%receive antennas, transmit antennas, taps, time axis
h = zeros(Mr,Nt,L,length(slow_time));
parfor tidx = 1:length(slow_time)
    temp = zeros(Mr,Nt,L);
    for lidx = 1:L
        var_lidx = Pt*1/tau_rms*ts*exp(-1/tau_rms*(lidx-1)*ts);
        D_fd = diag(exp(1j*(2*pi*fdmax*transpose(cosd(AoA))*(tidx-1)*TS+phi(:,lidx))));

        %         sqrt(sigma_lidx^2/N)* a_rx*D_fd*a_tx
        temp(1:Mr,1:Nt,lidx)  = sqrt(var_lidx/N)*a_rx*D_fd*a_tx;

    end
    h(:,:,:,tidx) = temp;
end
label_str = "h11";
semilogy(slow_time,squeeze(abs(h(1,1,1,:))),'-o',"DisplayName",label_str,'MarkerIndices',1:floor(length(slow_time)/5):length(slow_time),LineWidth=0.5)
hold on
label_str = "h12";
semilogy(slow_time,squeeze(abs(h(1,2,1,:))),'-x',"DisplayName",label_str,'MarkerIndices',1:floor(length(slow_time)/5):length(slow_time),LineWidth=0.5)
label_str = "h21";
semilogy(slow_time,squeeze(abs(h(2,1,1,:))),'-^',"DisplayName",label_str,'MarkerIndices',1:floor(length(slow_time)/5):length(slow_time),LineWidth=0.5)
label_str = "h22";
semilogy(slow_time,squeeze(abs(h(2,2,1,:))),'-*',"DisplayName",label_str,'MarkerIndices',1:floor(length(slow_time)/5):length(slow_time),LineWidth=0.5)
grid on
lgd = legend();
set(lgd,'Interpreter','latex');
set(lgd,'FontSize',12);
ylim([1e-3,100])
xlabel("time")
ylabel("amplitude")
title("Time Varying Frequency Selective 2x2 MIMO channel DP:"+string(tau_rms*1e6)+"us, fdmax:"+string(fdmax))

for k = 1:2
    for l = 1:2
        for i=1:2
            for j = 1:2
                c = corrcoef(squeeze(abs(h(i,j,1,:))),squeeze(abs(h(1,1,1,:))));
                c(1,2)
            end
        end
    end
end


h11 = squeeze(abs(h(1,1,1,:)));
h12 = squeeze(abs(h(1,2,1,:)));
h21 = squeeze(abs(h(2,1,1,:)));
h22 = squeeze(abs(h(2,2,1,:)));
% corrcoef(squeeze(abs(h(2,1,1,:))),squeeze(abs(h(2,2,1,:))))
% corrcoef(squeeze(abs(h(1,1,1,:))),squeeze(abs(h(2,2,1,:))))
%correlation between the antennas on the first tap
c = real(corrcoef([h11,h12,h21,h22]))

save("h.mat","h","fs","ts","T","t","fdmax","L","slow_time","TS","tau_rms")

%delay spread to frequency domain
figure;
%subcarriers
N = 128;
h_block = h(:,:,:,100);
h11 = squeeze(h_block(1,1,:));
h12 = squeeze(h_block(1,2,:));
h21 = squeeze(h_block(2,1,:));
h22 = squeeze(h_block(2,2,:));



H11 = fft(h11,N)/sqrt(N);
H12 = fft(h12,N)/sqrt(N);
H21 = fft(h21,N)/sqrt(N);
H22 = fft(h22,N)/sqrt(N);
%duration of the fast time axis, assuming N taps same as sub-carriers.
T = N*ts;
Delf = 1/T;
f_axis = -fs/2:Delf:fs/2-Delf;
label_str = "H11";
plot(1:N,10*log10(abs(H11)),'-o',"DisplayName",label_str,'MarkerIndices',1:10:N,LineWidth=0.5);
hold on
label_str = "H12";
plot(1:N,10*log10(abs(H12)),'-x',"DisplayName",label_str,'MarkerIndices',1:10:N,LineWidth=0.5);
label_str = "H21";
plot(1:N,10*log10(abs(H21)),'-^',"DisplayName",label_str,'MarkerIndices',1:10:N,LineWidth=0.5);
label_str = "H22";
plot(1:N,10*log10(abs(H22)),'-*',"DisplayName",label_str,'MarkerIndices',1:10:N,LineWidth=0.5);
xlabel("Sub-carrier number")
ylabel("Channel Amplitude")
grid on
lgd = legend();
set(lgd,'Interpreter','latex');
set(lgd,'FontSize',12);
title("Time Varying Frequency Selective 2x2 MIMO channel DP:"+string(tau_rms*1e6)+"us, fdmax:"+string(fdmax))
