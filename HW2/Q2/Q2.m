close all
clear all
clc




%mmse beamforming without interferers Npilots   

Ns = 100;
%assume NS/2 are pilots
Npilotsvec = 5:50;
M = 1;
N = 4;
c = 3e8;
f = 9e9;
lambda = c/f;
d = lambda/2;


AoA = 45;
AoD = 0;

ber = zeros(length(Npilotsvec),1);
ber2 = zeros(length(Npilotsvec),1);
ber3 = zeros(length(Npilotsvec),1);

Ntrials = 10000;
a0 = array_response(AoA,N,lambda,d);
parfor idy = 1:length(Npilotsvec)
    Npilots = Npilotsvec(idy);
    for n = 1:Ntrials
        %generate Npilots and (Ns-Npilots) data symbols 
        data = randi([0,3],[Ns,1]);
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
        snr_db = -1;

        
        %add noise
        [rx_frame,ns,txx] = add_noise(rx_frame,snr_db);
        
        rx_pilots = rx_frame(:,1:Npilots);
        rx_data = rx_frame(:,Npilots+1:end);
        
        
        Rrr = (rx_pilots)*(rx_pilots)'; 
        
        
        %idea MMSE weight
        w =pinv(Rrr)*a0;
        %mmse with channel estimation
        w2 = pinv(Rrr)*rx_pilots*conj(tx_frame(1:Npilots))/Npilots;
        %Max SNR beamformer
        w3 = a0;

        z_data = w'*rx_data;
        z_data2 = w2'*rx_data;
        z_data3 = w3'*rx_data;
        data_rx = qamdemod(z_data,4)';
        data_rx2 = qamdemod(z_data2,4)';
        data_rx3 = qamdemod(z_data3,4)';
        
        ber(idy) = ber(idy) +  1/Ntrials* (100-100*sum(data(Npilots+1:end)==data_rx)/(Ns-Npilots));
        ber2(idy) = ber2(idy) +  1/Ntrials* (100-100*sum(data(Npilots+1:end)==data_rx2)/(Ns-Npilots));
        ber3(idy) = ber3(idy) +  1/Ntrials* (100-100*sum(data(Npilots+1:end)==data_rx3)/(Ns-Npilots));

    end
end


% beampattern(w,N,lambda,d)
figure;
semilogy(Npilotsvec,ber,'-o')
hold on

grid on
xlabel('Npilots')
ylabel('BER');
semilogy(Npilotsvec,ber2,'-+')
semilogy(Npilotsvec,ber3,'-p')
legend("MMSE-perfectCSI","MMSE-estCSI","Max SNR")
title("Npilots - MMSE-perfectCSI, MMSE-estCSI, Max SNR ")


%%

%mmse beamforming without interferers SNR vec 

Ns = 100;
%assume NS/2 are pilots
Npilots = 50;
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
ber2 = zeros(length(SNR_vec),1);
ber3 = zeros(length(SNR_vec),1);

Ntrials = 10000;
a0 = array_response(AoA,N,lambda,d);
for idy = 1:length(SNR_vec)
    for n = 1:Ntrials
        %generate Npilots and (Ns-Npilots) data symbols 
        data = randi([0,3],[Ns,1]);
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
        
        rx_pilots = rx_frame(:,1:Npilots);
        rx_data = rx_frame(:,Npilots+1:end);
        
        
        Rrr = (rx_pilots)*(rx_pilots)'; 
        
        
        %idea MMSE weight
        w = pinv(Rrr)*a0;
        %mmse with channel estimation
        w2 = pinv(Rrr)*rx_pilots*conj(tx_frame(1:Npilots))/Npilots;
        %Max SNR beamformer
        w3 = a0;

        z_data = w'*rx_data;
        z_data2 = w2'*rx_data;
        z_data3 = w3'*rx_data;
        data_rx = qamdemod(z_data,4)';
        data_rx2 = qamdemod(z_data2,4)';
        data_rx3 = qamdemod(z_data3,4)';
        
        ber(idy) = ber(idy) +  1/Ntrials* (100-100*sum(data(Npilots+1:end)==data_rx)/(Ns-Npilots));
        ber2(idy) = ber2(idy) +  1/Ntrials* (100-100*sum(data(Npilots+1:end)==data_rx2)/(Ns-Npilots));
        ber3(idy) = ber3(idy) +  1/Ntrials* (100-100*sum(data(Npilots+1:end)==data_rx3)/(Ns-Npilots));

    end
end


beampattern(w2,N,lambda,d,true)
beampattern(w3,N,lambda,d,true)
figure;
% semilogy(SNR_vec,ber,'-o')
hold on

grid on
xlabel('SNR(dB)')
ylabel('BER');
semilogy(SNR_vec,ber2,'-+')
semilogy(SNR_vec,ber3,'-p')
legend("MMSE-estCSI","Max SNR")
title("SNR vs BER - MMSE-estCSI, Max SNR ")

%%


Ns = 1000;
%assume NS/2 are pilots
Npilots = 500;
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

%MMSE
ber_noint = zeros(length(SNR_vec),1);
ber_1int = zeros(length(SNR_vec),1);
ber_2int = zeros(length(SNR_vec),1);

ber_max_snr_noint = zeros(length(SNR_vec),1);
ber_max_snr_1int = zeros(length(SNR_vec),1);
ber_max_snr_2int = zeros(length(SNR_vec),1);

Ntrials = 10000;
%AoA for interferer
AoA_int = 0;
AoA_int2 = 25;
a0 = array_response(AoA,N,lambda,d);
a1 = array_response(AoA_int,N,lambda,d);
a2 = array_response(AoA_int2,N,lambda,d);
for idy = 1:length(SNR_vec)
    for n = 1:Ntrials
        %generate Npilots and (Ns-Npilots) data symbols 
        data = randi([0,3],[Ns,1]);
        iin = randi([0,3],[Ns,1]);
        iin2 = randi([0,3],[Ns,1]);
        iin_sym2 =  1/sqrt(2)*qammod(iin2,4);
        iin_sym =  1/sqrt(2)*qammod(iin,4);
        
        
        % modulate the data symbols
        syms = 1/sqrt(2)*qammod(data,4);
        
        
        %transmit symbols along AoD
        w = ones(M,1);
        tx_frame = transmit(syms,w,AoD,lambda,d);
        
        
         
        % NxNs receive matrix
        rx_frame_noint = a0*transpose(tx_frame);
        rx_frame_1int = rx_frame_noint+1/(sqrt(2))*a1*transpose(iin_sym);
        rx_frame_2int = rx_frame_noint+1/(sqrt(2))*a1*transpose(iin_sym)+1/(sqrt(2))*a2*transpose(iin_sym2);
        

        snr_db = SNR_vec(idy);
        %add noise
        [rx_frame_1int,ns,txx] = add_noise(rx_frame_1int,snr_db);
        rx_frame_noint = rx_frame_noint + ns;
        rx_frame_2int = rx_frame_2int + ns;


        rx_pilots_noint = rx_frame_noint(:,1:Npilots);
        rx_pilots_1int = rx_frame_1int(:,1:Npilots);
        rx_pilots_2int = rx_frame_2int(:,1:Npilots);
        rx_data_noint = rx_frame_noint(:,Npilots+1:end);
        rx_data_1int = rx_frame_1int(:,Npilots+1:end);
        rx_data_2int = rx_frame_2int(:,Npilots+1:end);
        
        Rrr_noint = rx_pilots_noint*rx_pilots_noint';
        Rrr_1int = (rx_pilots_1int)*(rx_pilots_1int)'; 
        Rrr_2int = (rx_pilots_2int)*(rx_pilots_2int)'; 
        
        %MMSE weight with CSI estimation
        w_noint =  pinv(Rrr_noint)*rx_pilots_noint*conj(tx_frame(1:Npilots))/Npilots;
        w_1int =  pinv(Rrr_1int)*rx_pilots_1int*conj(tx_frame(1:Npilots))/Npilots;
        w_2int =  pinv(Rrr_2int)*rx_pilots_2int*conj(tx_frame(1:Npilots))/Npilots;
        %Maxsnr
        w_maxsnr = a0;
        
        z_data_noint = w_noint'*(rx_data_noint);
        z_data_1int = w_1int'*rx_data_1int;
        z_data_2int = w_2int'*rx_data_2int;
        z_data_maxsnr_noint = w_maxsnr'*rx_data_noint;
        z_data_maxsnr_1int = w_maxsnr'*rx_data_1int;
        z_data_maxsnr_2int = w_maxsnr'*rx_data_2int;

        data_rx_noint =  qamdemod(z_data_noint,4)';
        data_rx_1int = qamdemod(z_data_1int,4)';
        data_rx_2int = qamdemod(z_data_2int,4)';
        
        data_rx_max_snr_noint = qamdemod(z_data_maxsnr_noint,4)'; 
        data_rx_max_snr_1int = qamdemod(z_data_maxsnr_1int,4)'; 
        data_rx_max_snr_2int = qamdemod(z_data_maxsnr_2int,4)'; 


        ber_noint(idy) = ber_noint(idy) +  1/Ntrials* (100-100*sum(data(Npilots+1:end)==data_rx_noint)/(Ns-Npilots));
        ber_1int(idy) = ber_1int(idy) +  1/Ntrials* (100-100*sum(data(Npilots+1:end)==data_rx_1int)/(Ns-Npilots));
        ber_2int(idy) = ber_2int(idy) +  1/Ntrials* (100-100*sum(data(Npilots+1:end)==data_rx_2int)/(Ns-Npilots));

        ber_max_snr_noint(idy) = ber_max_snr_noint(idy) +  1/Ntrials* (100-100*sum(data(Npilots+1:end)==data_rx_max_snr_noint)/(Ns-Npilots));
        ber_max_snr_1int(idy) = ber_max_snr_1int(idy) +  1/Ntrials* (100-100*sum(data(Npilots+1:end)==data_rx_max_snr_1int)/(Ns-Npilots));
        ber_max_snr_2int(idy) = ber_max_snr_2int(idy) +  1/Ntrials* (100-100*sum(data(Npilots+1:end)==data_rx_max_snr_2int)/(Ns-Npilots));


    end
end

beampattern(w_noint,N,lambda,d,true)
beampattern(w_1int,N,lambda,d,true)
beampattern(w_2int,N,lambda,d,true)
beampattern(a0,N,lambda,d,true)

% beampattern(w,N,lambda,d,true);
figure;

plot(SNR_vec,ber_noint,'-o')
hold on
plot(SNR_vec,ber_1int,'-o')
plot(SNR_vec,ber_2int,'-o')

plot(SNR_vec,ber_max_snr_noint,'-x')
plot(SNR_vec,ber_max_snr_1int,'-x')
plot(SNR_vec,ber_max_snr_2int,'-x')
title("BER: MMSE vs Max SNR for no int, 1 int, 2 int")
grid on
xlabel('SNR(dB)');
ylabel('BER');
legend("MMSE - no int","MMSE - 1 int","MMSE - 2 int","MaxSNR - no int","MaxSNR - 1 int","MaxSNR - 2 int")