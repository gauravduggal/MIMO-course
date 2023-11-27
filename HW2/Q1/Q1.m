close all
clear all
clc


    

Ns = 100;
%assume NS/2 are pilots
Npilots = floor(Ns/2);
M = 1;
N = 4;
c = 3e8;
f = 9e9;
lambda = c/f;
d = lambda/2;


AoA = 0;
AoD = 0;
%SNR per receive antenna
SNR_vec = 5:-0.5:-1;
ber = zeros(length(SNR_vec),1);
Ntrials = 1000;

a0 = array_response(AoA,N,lambda,d);
for idy = 1:length(SNR_vec)
    for n = 1:Ntrials
        %generate Npilots and (Ns-Npilots) data symbols 
        data = randi([0,3],[Ns,1]);
 

        data(1:Npilots) = 0;
        % temp =repelem([0;3],1,Npilots);
        % temp = temp(:);
        % data(1:Npilots) = temp(1:Npilots);
        
        
        % modulate the data symbols
        syms = 1/sqrt(2)*qammod(data,4);
        
        
        %transmit symbols along AoD
        w = ones(M,1);
        tx_frame = transmit(syms,w,AoD,lambda,d);
        
        
         
        % NxNs receive matrix
        rx_frame = a0*transpose(tx_frame);
        snr_db = SNR_vec(idy);

        
        %add noise
        [rx_frame,ns,txx] = add_noise(rx_frame,snr_db);
        
      
        
        %max SNR weight
        w = a0;
        
        z_data = w'*rx_frame;
        
        data_rx = qamdemod(z_data,4)';
        ber(idy) = ber(idy) +  1/Ntrials* (100-100*sum(data==data_rx)/(Ns));
    end
end


% beampattern(w,N,lambda,d)
figure;
semilogy(SNR_vec,ber)
grid on
xlabel('SNR(dB)')
ylabel('BER');
hold on

% bit error rate for 1 tx- 1 rx for QPSK assuming SNR is 
ber_1tx_1rx = 100*qfunc(sqrt(10.^(SNR_vec./10)));
semilogy(SNR_vec,ber_1tx_1rx,'-o')
xlabel("SNR (dB)")
ylabel("BER %")
title("BER for TX-RX beamforming")

grid on
ylim([1e-2 20])




%%




    

Ns = 1000;

M = 4;
N = 1;
c = 3e8;
f = 9e9;
lambda = c/f;
d = lambda/2;


AoA = 0;
AoD = 0;
%SNR per receive antenna
SNR_vec = 5:-0.5:-1;
ber = zeros(length(SNR_vec),1);
Ntrials = 100;
% rx_frame = zeros(N,1);
a0 = array_response(AoA,M,lambda,d);
for idy = 1:length(SNR_vec)
    for n = 1:Ntrials
        %generate Npilots and (Ns-Npilots) data symbols 
        data = randi([0,3],[Ns,1]);
 

        
        
        % modulate the data symbols
        syms = 1/sqrt(2)*qammod(data,4);
        
        
        %transmit symbols along AoD
        w = ones(M,1);
        tx_frame = a0*transpose(transmit(syms,w,AoD,lambda,d));
        
        
         
        % NxNs receive matrix
        rx_frame = tx_frame/M;
        snr_db = SNR_vec(idy);

        
        %add noise
        [rx_frame,ns,txx] = add_noise(rx_frame,snr_db);
        

        
           
        
        %max SNR weight
        w = a0;
        
        z_data = w'*rx_frame;
        
        data_rx = qamdemod(z_data,4)';
        ber(idy) = ber(idy) +  1/Ntrials* (100-100*sum(data==data_rx)/(Ns));
    end
end


% beampattern(w,N,lambda,d)
% figure;
semilogy(SNR_vec,ber,'square')
grid on
% title("BER for TX beamforming")
xlabel('SNR(dB)')
ylabel('BER');
legend("1tx-4rx","1tx-1rx theoretical","4tx-1rx")

%% Q1(c)

M = 4;
N = 1;
c = 3e8;
f = 9e9;
lambda = c/f;
d = lambda/2;
AoD = 35;
theta_vec = -90:1:90;
af = abs(array_factor(theta_vec-AoD, d, lambda,M));
figure;

plot(theta_vec,af);
xlabel("AoD wrt to boresight (degrees)->")
ylabel("Array Gain")
title("TX side beamforming at A0D 35 degrees")


%% Q1 (d). 
M = 1;
N = 4;
c = 3e8;
f = 9e9;
lambda = c/f;
d = lambda/2;
AoA = 15;
theta_vec = -90:1:90;

af1 = abs(array_factor(theta_vec-AoA, d, lambda,N));
figure;

plot(theta_vec,af1);
hold on

AoA = 65;
theta_vec = -90:1:90;
af2 = abs(array_factor(theta_vec-AoA, d, lambda,N));
plot(theta_vec,af2);
legend("AoA 15 degrees","AoA 65 degrees")
xlabel("AoD wrt to boresight (degrees)->")
ylabel("Array Gain")
title("RX beamsteering for AOA 15 degrees and 65 degrees")