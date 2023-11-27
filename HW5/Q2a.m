clear all
clc
close all


Ntrials = 1000;

snrvec = 8:-1:0;
ber =  zeros(size(snrvec));
%16-QAM
k = 4;
% sumber of subcarriers
N = 128;
Dn = dftmtx(N)/sqrt(N);
for snridx=1:length(snrvec)
    tic
    snr = snrvec(snridx);
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
%         tx = (Dn')*sym;

%         sum(abs(tx1-tx))
        
        rx = add_noise_td(tx,snr,k);
        
        %magic sauce section
        s_est = fft(rx)/sqrt(N);
%         s_est = Dn*rx;
        s_est_dec = qamdemod(s_est,2^k,'UnitAveragePower',true);

        ber_snr = ber_snr +  1/Ntrials*(1 - 1/(k*N)*sum(sum(int2bit(s_est_dec,k)==int2bit(sym_dec,k))));

    end
    ber(snridx) = ber_snr;
    toc
end
legend_str = "OFDM-16QAM-AWGN";
semilogy(snrvec,ber,"-x","DisplayName",legend_str,LineWidth=1.5);
ylim([1e-3 1])
xlim([0,8])
hold on
legend


parfor snridx=1:length(snrvec)
    tic
    snr = snrvec(snridx);
    ber_snr = 0;
    for trials=1:Ntrials
        Nbits = k;
        data_bits = floor(2*rand(Nbits,1));
%         bits = data_bits;
        temp = reshape(data_bits,[k,1]);
        sym_dec = bit2int(temp,k);
        %assume symbols are in frequency domain
        sym = qammod(sym_dec,2^k,'UnitAveragePower',true);
        %modulate the symbols and convert to time domain of N samples
        
%         rx = sym;
        rx = add_noise(sym,snr,k);
        
        %magic sauce section
        s_est = rx;
        s_est_dec = qamdemod(s_est,2^k,'UnitAveragePower',true);

        ber_snr = ber_snr +  1/Ntrials*(1 - 1/(k)*sum(sum(int2bit(s_est_dec,k)==int2bit(sym_dec,k))));

    end
    ber(snridx) = ber_snr;
    toc
end
legend_str = "Theory-16-QAM-AWGN";
semilogy(snrvec,ber,"-o","DisplayName",legend_str, LineWidth=1.5);
grid on
xlabel("E_b/N_0 (dB)")
ylabel("BER")
legend
