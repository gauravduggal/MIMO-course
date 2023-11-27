function [rx_sym] = add_noise_td2(tx_samples,snr_db,N_sc,ofdm_symbols)
snr = 10^(snr_db/10);
%number of time samples
N = size(tx_samples,1);
noise = randn(N,1) + 1j*randn(N,1);
%each ofdm sylmbol has N time samples where N is the number of sub carriers
%here snr is per subcarrier and each subcarrier has a symbol of unit power
%hence the normalisation is just with snr. The 2 is for I and Q samples
noise = noise/sqrt(2*snr);
% N
% sum(abs(tx_samples).^2)
% sum(abs(noise).^2)
% [10*log10(sum(abs(tx_samples).^2)/sum(abs(noise).^2)), snr_db]
%  snr_db
rx_sym = tx_samples+noise;

end

