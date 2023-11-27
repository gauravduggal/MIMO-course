close all
clear all
clc

N = 128;
Ntrials = 1000;
tau_rms = 100e-9;
L = 10;
fs = 100e6;
snrvec =35:-1:5;
ber = zeros(size(snrvec));
cp_n = 11;
k = 6;
H_trials = zeros(N,Ntrials);
for snridx = 1:length(snrvec)
    snr_db = snrvec(snridx);
    for trial = 1:Ntrials
        N = 128;

        
        Nbits = N*k;

        h = transpose(get_freq_selective_channel(tau_rms,fs,L));
        h = h/sum(abs(h).^2);
        H = fft(h,N)/sqrt(N);
        H_trials(:,trial) = sqrt(N*10^(snr_db/10))*H;
        dataIn = floor(2*rand(Nbits,1));

        syms = reshape(dataIn,[k,N]);
        sym_dec = bit2int(syms,k)';
        sym_tx = qammod(sym_dec,2^k,'UnitAveragePower',true);

        tx = ifft(sym_tx,N)*sqrt(N);

        cp = tx(N-cp_n+1:N,1);
        tx_cp = 1/sqrt((cp_n+N)/N)*[cp;tx];
%         tx_cp = [cp;tx];
        % plot(abs(tx))
        rx = conv(tx_cp,h);
        rx = add_noise_td2(rx,snr_db,1,1);
        % rx = rx;

        rx_sym = fft(rx(cp_n+1:N+cp_n,1),N)/sqrt(N);
        %%
        rx_sym = rx_sym./(sqrt(N)*H);

        rx_demod = qamdemod(rx_sym,2^k,'UnitAveragePower',true);

        dataOut  = int2bit(rx_demod,k);


        ber(snridx) = ber(snridx) + 1/Ntrials* sum(dataOut~=dataIn)/(N*k);
    end
end
semilogy(snrvec,ber)
hold on
semilogy(snrvec,get_ber_rayleigh(snrvec,k))
ylim([1e-4,1])
grid on
figure
histogram(10*log10(abs(H_trials(24,:).^2)))
grid on