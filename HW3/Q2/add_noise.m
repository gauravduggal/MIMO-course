
function [rx_sym] = add_noise(tx_sym,snr_db,k)
snr = 10^(snr_db/10);
L = size(tx_sym,1);
noise = randn(L,1) + 1j*randn(L,1);
noise = noise / sqrt(2*k*snr);
rx_sym = tx_sym+noise;

end