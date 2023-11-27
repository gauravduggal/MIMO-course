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


AoA = 45;
AoD = 0;
%SNR per receive antenna
SNR_vec = 5:-0.5:-1;
ber = zeros(length(SNR_vec),1);
ber_noint = zeros(length(SNR_vec),1);
Ntrials = 10000;
%AoA for interferer
AoA_int = 25;
a0 = array_response(AoA,N,lambda,10*lambda);
a1 = array_response(AoA_int,N,lambda,10*lambda);
for idy = 1:length(SNR_vec)
    for n = 1:Ntrials
        %generate Npilots and (Ns-Npilots) data symbols 
        data = randi([0,3],[Ns,1]);
        iin = randi([0,3],[Ns,1]);
        iin_sym =  1/sqrt(2)*qammod(iin,4);
        data(1:Npilots) = 0;

        
        % modulate the data symbols
        syms = 1/sqrt(2)*qammod(data,4);
        
        
        %transmit symbols along AoD
        w = ones(M,1);
        tx_frame = transmit(syms,w,AoD,lambda,d);
        
        
         
        % NxNs receive matrix
        rx_frame_noint = a0*transpose(tx_frame);
        rx_frame = rx_frame_noint+a1*transpose(iin_sym);
        snr_db = SNR_vec(idy);

        
        %add noise
        [rx_frame,ns,txx] = add_noise(rx_frame,snr_db);
        rx_frame_noint = rx_frame_noint + ns;
        %add rayleigh fading
        h = 1/sqrt(2)*(randn(N,1)+1j*randn(N,1));
        rx_frame = h.*rx_frame;
        rx_frame_noint = h.*rx_frame_noint;

        rx_pilots_noint = rx_frame_noint(:,1:Npilots);
        rx_pilots = rx_frame(:,1:Npilots);
        rx_data_noint = rx_frame_noint(:,Npilots+1:end);
        rx_data = rx_frame(:,Npilots+1:end);
        
        Rrr_noint = rx_pilots_noint*rx_pilots_noint';
        Rrr = (rx_pilots)*(rx_pilots)'; 
        
        %MMSE weight
        w_noint =  pinv(Rrr_noint)*rx_pilots_noint*conj(tx_frame(1:Npilots))/Npilots;
        w = pinv(Rrr)*rx_pilots*conj(tx_frame(1:Npilots))/Npilots;
%         w = a0;
        
        z_data_noint = w_noint'*(rx_data_noint);
        z_data = w'*rx_data;
        
        data_rx_noint =  qamdemod(z_data_noint,4)';
        data_rx = qamdemod(z_data,4)';
        ber(idy) = ber(idy) +  1/Ntrials* (100-100*sum(data(Npilots+1:end)==data_rx)/(Ns-Npilots));
        ber_noint(idy) = ber_noint(idy) +  1/Ntrials* (100-100*sum(data(Npilots+1:end)==data_rx_noint)/(Ns-Npilots));
    end
end


% beampattern(w,N,lambda,d,true);
figure;

plot(SNR_vec,ber_noint)

grid on
xlabel('SNR(dB)');
ylabel('BER');
% title("BER vs SNR for rayleigh with 1 interferer")

% ylim([30,40])
hold on
% figure;
plot(SNR_vec,ber)
xlabel('SNR(dB)');
ylabel('BER %');
title("BER vs SNR for rayleigh fading")
grid on
legend("no int","1 int")
% ylim([1e-2,10^2])