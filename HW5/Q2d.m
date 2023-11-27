clear all
clc
close all

%sampling frequency
fs = 100e6;
ts=1/fs;
Ntrials = 1000;
%delay spread
tau_rms = 100e-9;
%cycling prefix vector
cpvec = [4,8,16,32,64,127];
snrvec = 20:-0.5:0;
ber =  zeros(size(snrvec));
%BPSK
k = 1;
% sumber of subcarriers
N = 128;
Dn = dftmtx(N)/sqrt(N);
for cpidx = 1:length(cpvec)
    cp_n = cpvec(cpidx);
    parfor snridx=1:length(snrvec)
        tic
        snr = snrvec(snridx);
        ber_snr = 0;
        for trials=1:Ntrials
            %create 2 packets of bits
            Nbits = 2*N*k;
            data_bits_p1 = floor(2*rand(Nbits/2,1));
            data_bits_p2 = floor(2*rand(Nbits/2,1));

            temp = reshape(data_bits_p1,[k,N]);
            sym_dec_p1 = bit2int(temp,k)';
            temp = reshape(data_bits_p2,[k,N]);
            sym_dec_p2 = bit2int(temp,k)';

            %assume symbols are in frequency domain
            sym_p1 = qammod(sym_dec_p1,2^k,'UnitAveragePower',true);
            sym_p2 = qammod(sym_dec_p2,2^k,'UnitAveragePower',true);

            %modulate the symbols across subcarriers and convert to time domain with N time samples
            tx_p1 = conj(Dn)*sym_p1;
            tx_p2 = conj(Dn)*sym_p2;

            %create channel in time domain
            h = get_freq_selective_channel(tau_rms,fs,N);
            h = transpose(h);
            %create CP
            cp1 = tx_p1(N-cp_n:N,1);
            cp2 = tx_p2(N-cp_n:N,1);
            %append tx with cp
            tx = [cp1;tx_p1;cp2;tx_p2];

            rx_cp = conv(tx,h);

            %         rx_cp = add_noise_td(rx_cp,snr,k);
            %remove cycling prefix and tail from second packet
            
            rx = rx_cp(cp_n+N+1+cp_n+1+1:cp_n+N+1+cp_n+1+N);
            rx = add_noise_td(rx,snr,k);
            %Equalization
            %frequency domain channel
%             H = (fft(h)/sqrt(N));
            H = Dn*h;
            s_est = (Dn*rx)./(H);
            %         s_est = fft(rx);
            s_est_dec = qamdemod(s_est,2^k,'UnitAveragePower',true);

            ber_snr = ber_snr +  1/Ntrials*(1 - 1/(k*N)*sum(sum(int2bit(s_est_dec,k)==int2bit(sym_dec_p2,k))));

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
title("OFDM frequency selective fading")
theory_ber = 0.25./(10.^(snrvec/10));
legend_str = "OFDM-theory-flat-fading";
semilogy(snrvec,theory_ber,"-x","DisplayName",legend_str,LineWidth=1.5);
legend(Location="best")
xlabel("E_b/N_0 (dB)")
ylabel("BER")
ylim([1e-3,1])
% grid on
