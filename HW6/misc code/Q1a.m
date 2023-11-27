clear all
clc
close all


%sampling frequency
fs = 100e6;
ts=1/fs;
Ntrials = 100;
%delay spread
tau_rms = 100e-9;
%number of pilots
npvec = [16,32,64];
snrvec = 20:-0.5:0;
ber =  zeros(size(snrvec));
mse_channel_estimation = zeros(length(snrvec),length(npvec));
%QPSK
k = 1;
% sumber of subcarriers
N = 128;
Dn = dftmtx(N)/sqrt(N);
%cyclic prefix > delay spread
cp_n = 15;
%no of ofdm symbols in one packet assuming block fading i.e. prop to
%coherence time
ofdm_symbols = 1;
%no of channel taps
L = 10;


for npidx = 1:length(npvec) 
%pilot related variables 
Np = npvec(npidx);
pilot_spacing = floor(N/Np);
%pilot symbols to subcarrier indices
pilots_symbols_idx = k+1:pilot_spacing:N-pilot_spacing;
Np = length(pilots_symbols_idx);
%dft matrices
QNp = dftmtx(Np)/sqrt(Np);
QNpL = QNp(:,1:Np);
%data symbols to subcarrier indices
data_symbols_idx = 1:N;
data_symbols_idx(pilots_symbols_idx) = [];
Ns = length(data_symbols_idx);
%total bits over all packets including both pilots and data symbols
Nbits = ofdm_symbols*N*k;
temp = reshape(1:Nbits,[N*k,ofdm_symbols]);
pilot_bits = [];
data_bits=[];
for scidx=1:N-1
    %if subcarrier has a pilot symbol find bit indexes
    if sum(scidx==pilots_symbols_idx)
        if scidx==1
            bidx = 1:k;
        elseif scidx == N
            bidx = (N-1)*k+1:N*k;
        else
            bidx = scidx*k:(scidx+1)*k-1; 
        end
        pilot_bits = [pilot_bits, bidx];
    else
        if scidx==1
            bidx = 1:k;
        elseif scidx == N
            bidx = (N-1)*k+1:N*k;
        else
            bidx = scidx*k:(scidx+1)*k-1; 
        end
        data_bits = [data_bits, bidx];
    end
%     deletions = [ deletions , pilots_symbols_idx(pidx)*k:pilots_symbols_idx(pidx)*k+k-1];
end
temp2 = temp;
temp2(data_bits)=[];
temp(pilot_bits,:)=[];
data_bits_idx = temp(:);
pilot_bits_idx = temp2(:);
%%

parfor snridx=1:length(snrvec)
    tic
    snr = snrvec(snridx);
    ber_snr = 0;
    mse_channel_estimation_snr = 0;
    tx_pilots = zeros(ofdm_symbols,Ns);
    for trials=1:Ntrials
        mse_channel_estimation_trial = 0;
        dataEnc  = floor(2*rand(Nbits,1));
        temp = reshape(dataEnc,[k,ofdm_symbols*N]);
        sym_dec = transpose(bit2int(temp,k));

            
        %modulate symbols across subcarriers,
        sym = qammod(sym_dec,2^k,'UnitAveragePower',true);
        %add cyclic prefix

        %assume symbols are in frequency domain and convert to time domain with N time samples
        tx = zeros(N*ofdm_symbols,1);
        for sidx = 1:ofdm_symbols
            st_idx = (sidx-1)*N+1;
            en_idx = st_idx+N-1;
            tx(st_idx:en_idx) = ifft(sym(st_idx:en_idx,1))*sqrt(N);
            
            tx_pilots(sidx,1:Np) = sym(pilots_symbols_idx+(sidx-1)*N);
        end
        %yucky matlab syntax section to add cyclic prefix
        temp = reshape(tx,[N,ofdm_symbols]);
        temp = transpose(temp);

        cp = temp(:,N-cp_n+1:N);
        temp2 = [cp,temp];
        %dont ask why is there a transpose, stupid matlab syntax
        tx_cp = reshape(transpose(temp2),[(N+cp_n)*ofdm_symbols,1]);

        h = get_freq_selective_channel(tau_rms,fs,Np);
        h = transpose(h);

        %this results in packets*(N+CP_N)+N-1 samples, the last N-1 samples is the
        %tail
        rx_cp = conv(tx_cp,h);
        rx_cp = add_noise_td(rx_cp,snr,k);

        %frequency domain channel
        H = fft(h,N)/sqrt(N);


        temp = zeros(ofdm_symbols,N);
        encoded_received_bits = zeros(ofdm_symbols*Ns*k,1);
        for sidx = 1:ofdm_symbols
            %strip cyclic prefix from every packet and the tail from 12th
            %packet to get integration length for each OFDM symbol
            st_idx = (sidx-1)*N+(sidx)*cp_n+1;
            en_idx = st_idx+N-1;
            %get integration length
            temp(sidx,:) = rx_cp(st_idx:en_idx,1);
            %matched filter  to convert to received distorted symbols in freq domain
            temp(sidx,:) = fft(temp(sidx,:),N)/(sqrt(N));
            
            %channel estimation part
            %pick pilots in ofdm symbol one
            ep = temp(sidx,pilots_symbols_idx);
            %LS channel estimation
            h_est = (QNpL'*QNpL)^-1*QNpL'*inv(diag(tx_pilots(sidx,1:Np)))*transpose(ep)/sqrt(Np);
            %freq domain of estimated channel
            H_est = fft(h_est,N)/sqrt(N);

            
            %pick received symbols (pilots + data)
            received_data_symbols = temp(sidx,:);
            %equalisation with estimated channel
            received_data_symbols = inv(diag(H_est))*transpose(received_data_symbols);
            %demodulate using modulation scheme
            s_est_dec = qamdemod(received_data_symbols,2^k,'UnitAveragePower',true);
            %pick data symbols
            s_est_dec = s_est_dec(data_symbols_idx);
            packet_idx = (sidx-1)*k*Ns+1:(sidx)*k*Ns;
            encoded_received_bits(packet_idx,1) = int2bit(s_est_dec,k);
            mse_channel_estimation_trial = mse_channel_estimation_trial + 1/ofdm_symbols*sum(abs(h-h_est).^2);

        end
        %             bit_error_current_packet = 1 - 1/(k*N)* sum(sum(int2bit(s_est_dec,k)==int2bit(sym_dec(packet_idx,1),k)));
        %             ber_snr = ber_snr +  1/(packets*Ntrials)*bit_error_current_packet;
        dataOut = encoded_received_bits;
        ber_snr = ber_snr + 1/Ntrials * (1-(sum((dataEnc(data_bits_idx)==dataOut))/(Ns*k*ofdm_symbols)));
        mse_channel_estimation_snr = mse_channel_estimation_snr + 1/Ntrials *mse_channel_estimation_trial;
        %             ber_snr = ber_snr + 1/Ntrials*(1 - sum(dataOut == dataIn)/Nbits);
        %                              ber_snr = ber_snr + 1/Ntrials*(1 - sum(dataEnc == encoded_received_bits)/3072);
        %             ber_snr = ber_snr + 1/Ntrials*(errorStats(1));
    end
    ber(snridx) = ber_snr;
    mse_channel_estimation(snridx,npidx) = mse_channel_estimation_snr;
    toc
end
legend_str ="Delay spread:"+string(floor(tau_rms*1e9))+"ns, Cyclic prefix:"+string(floor(cp_n*ts*1e9))+"ns,"+" Pilots:"+string(npvec(npidx));
semilogy(snrvec,ber,"-x","DisplayName",legend_str,LineWidth=1.5);

hold on
grid on
title("OFDM-BER in AWGN vs Npilots")
xlabel("E_b/N_0 (dB)")
ylabel("BER")
ylim([1e-3,1])
end
% theory_ber = 0.25./(10.^(snrvec/10));
% legend_str = "OFDM-BER in AWGN vs Npilots";
% semilogy(snrvec,theory_ber,"-x","DisplayName",legend_str,LineWidth=1.5);
lgd = legend(Location="best");
set(lgd,'Interpreter','latex');
set(lgd,'FontSize',12);
grid on

figure;

for npidx = 1:length(npvec)
legend_str = "Pilots Np:"+string(npvec(npidx));
semilogy(snrvec,mse_channel_estimation(:,npidx),"-x","DisplayName",legend_str,LineWidth=1.5)
hold on
ylim([1e-2 10])
xlim([0,20])
end
grid on
ylabel("Error")
xlabel("E_b/N_0 (dB)")
lgd = legend(Location="best");
set(lgd,'Interpreter','latex');
set(lgd,'FontSize',12);
grid on
title("Mean Squared Error- Channel Estimation")

