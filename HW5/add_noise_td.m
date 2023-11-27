
function [rx_sym] = add_noise_td(tx_samples,snr_db,k)
snr = 10^(snr_db/10);
%number of time samples
N = size(tx_samples,1);
noise = randn(N,1) + 1j*randn(N,1);
noise = noise / sqrt(2*k*snr);
rx_sym = tx_samples+noise;

end