clear all
clc
close all





%sampling frequency
fs = 100e6;
ts=1/fs;
Ntrials = 300;
%delay spread
tau_rms = 100e-9;
%cycling prefix vector
dsvec = [10,100,1000]*1e-9;
snrvec = 20:-0.5:0;
ber =  zeros(size(snrvec));
%QPSK
k = 2;
% sumber of subcarriers
N = 128;
Dn = dftmtx(N)/sqrt(N);
for dsidx = 1:length(dsvec)
    tau_rms = dsvec(dsidx);
    %cyclic prefix > delay spread
    cp_n = 50;
    parfor snridx=1:length(snrvec)
        tic
        snr = snrvec(snridx);
        ber_snr = 0;
        for trials=1:Ntrials
                            %create channel in time domain
                %             Nbits = 12*floor(2/3*N*k);
                % bchnumerr(2047)
                % we get 2047 encoded bits, on running bchnumerr(2047) we
                % pick one value of K i.e. in this case 1409 which leads to
                % error correction capability of 60 bits
                Nbits = 1409;
                dataIn = floor(2*rand(Nbits,1));
                

                dataEnc = transpose(BCH_encode(dataIn,11,Nbits));
                %add one bit to make multiple of 8
                dataEnc  = [0;dataEnc];

                temp = reshape(dataEnc,[k,8*N]);
                sym_dec = transpose(bit2int(temp,k));

                %modulate symbols across subcarriers,
                sym = qammod(sym_dec,2^k,'UnitAveragePower',true);
                %add cyclic prefix

                %assume symbols are in frequency domain and convert to time domain with N time samples
                tx = [];
                for sidx = 1:8
                    st_idx = (sidx-1)*N+1;
                    en_idx = st_idx+N-1;
                    tx = [tx;ifft(sym(st_idx:en_idx,1))*sqrt(N)];
                    %                 tx = [tx;conj(Dn)*sym(st_idx:en_idx,1)];
                end
                %yucky matlab syntax section to add cyclic prefix
                temp = reshape(tx,[N,8]);
                temp = transpose(temp);

                cp = temp(:,N-cp_n+1:N);
                temp2 = [cp,temp];
                %dont ask why is there a transpose, stupid matlab syntax
                tx_cp = reshape(transpose(temp2),[(N+cp_n)*8,1]);

                H = (randn(1)+1j*randn(1))/(sqrt(2));
                
                %this results in 8*(N+CP_N)+N-1 samples, the last N-1 samples is the
                %tail
%                 rx_cp = conv(tx_cp,h);
                
                rx_cp = add_noise_td(H*tx_cp,snr,k);
                
                

                temp = zeros(8,N);
                encoded_received_bits = zeros(8*N*k,1);
                for sidx = 1:8
                    %strip cyclic prefix from every packet and the tail from 12th
                    %packet to get integration length for each OFDM symbol
                    st_idx = (sidx-1)*N+(sidx)*cp_n+1;
                    en_idx = st_idx+N-1;
                    temp(sidx,:) = rx_cp(st_idx:en_idx,1);
                    %demodulate symbols off each subcarrier
                    temp(sidx,:) = fft(temp(sidx,:))/(sqrt(N));
                    %                 temp(sidx,:) = Dn*temp(sidx,:);
                    %equalise channel - thanks to cyclic prefix the convolution
                    %results in an efective circular convolution therefore in
                    %the frequency domain we can just divide the freq domain
                    % tx signal by freq domain channelby the channel
                    temp(sidx,:) = temp(sidx,:)/(transpose(H));
                    %decode symbols
                    s_est_dec = qamdemod(temp(sidx,:),2^k,'UnitAveragePower',true)';
                    packet_idx = (sidx-1)*k*N+1:(sidx)*k*N;
                    packet_idx_rev = (sidx)*k*N:-1:(sidx-1)*k*N+1;
                    encoded_received_bits(packet_idx,1) = int2bit(s_est_dec,k);

                end
                %             bit_error_current_packet = 1 - 1/(k*N)* sum(sum(int2bit(s_est_dec,k)==int2bit(sym_dec(packet_idx,1),k)));
                %             ber_snr = ber_snr +  1/(8*Ntrials)*bit_error_current_packet;
                dataOut = transpose(BCH_decode(encoded_received_bits(2:end),11,Nbits));
                ber_snr = ber_snr + 1/Ntrials * (1-sum((dataIn==dataOut))/Nbits);
                %             ber_snr = ber_snr + 1/Ntrials*(1 - sum(dataOut == dataIn)/Nbits);
%                              ber_snr = ber_snr + 1/Ntrials*(1 - sum(dataEnc == encoded_received_bits)/3072);
%             ber_snr = ber_snr + 1/Ntrials*(errorStats(1));
        end
        ber(snridx) = ber_snr;
        toc
    end
    legend_str ="Delay spread:"+string(floor(tau_rms*1e9))+"ns, Cyclic prefix:"+string(floor(cp_n*ts*1e9))+"ns";
    semilogy(snrvec,ber,"-x","DisplayName",legend_str,LineWidth=1.5);
    % ylim([1e-5 1])
    % xlim([0,20])
    hold on
    grid on
end
title("OFDM flat fading with coding")
% theory_ber = 0.25./(10.^(snrvec/10));
% legend_str = "OFDM-theory-flat-fading";
% semilogy(snrvec,theory_ber,"-x","DisplayName",legend_str,LineWidth=1.5);
legend(Location="best")
xlabel("E_b/N_0 (dB)")
ylabel("BER")
ylim([1e-3,1])
% grid on
