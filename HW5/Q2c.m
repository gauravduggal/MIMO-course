clear all
clc
close all


Ntrials = 1000;

snr_db = 6;
snr = 10^(snr_db/10);
%BPSK
k = 1;
% sumber of subcarriers
N = 128;
fvec = 0:0.005:0.2-0.005;
ber =  zeros(size(fvec));

parfor fidx=1:length(fvec)
    tic
    f_offset = fvec(fidx);
    ber_snr = 0;
    for trials=1:Ntrials
        Nbits = N*k;
        data_bits = floor(2*rand(Nbits,1));
        bits = data_bits;
        temp = reshape(bits,[k,N]);
        sym_dec = bit2int(temp,k)';
        %assume symbols are in frequency domain
        sym = qammod(sym_dec,2^k,'UnitAveragePower',true);
        %modulate the symbols and convert to time domain of N samples
        tx = ifft(sym)*sqrt(N);
        
%         rx = tx;
        %time domain flat fading channel

        
        %apply channel in frequency domain
%         rx = ifft((fft(tx)/(sqrt(N))).*H);

        rx = add_noise_td(tx,snr,k);
        %apply frequency offset in time domain
        rx = rx.*exp(1j*2*pi*f_offset*(1:N)'/N);
        %magic sauce section
        s_est = (fft(rx)/sqrt(N));%./H;
        s_est_dec = qamdemod(s_est,2^k,'UnitAveragePower',true);

        ber_snr = ber_snr +  1/Ntrials*(1 - 1/(Nbits)*sum(sum(int2bit(s_est_dec,k)==int2bit(sym_dec,k))));

    end
    ber(fidx) = ber_snr;
    toc
end
legend_str = "OFDM-AWGN-BPSK-E_b/N_0=6dB";
semilogy(fvec,ber,"-x","DisplayName",legend_str,LineWidth=1.5);
ylim([1e-3 1])
hold on
legend
title("BER vs frequency offset")
xlabel("frequency offset (normalised frequency)")
ylabel("BER")
grid on


