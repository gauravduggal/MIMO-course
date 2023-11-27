
function [rx_sym] = add_noise(tx_sym,snr_db,k)
snr = 10^(snr_db/10);
Nr = size(tx_sym,1);
noise = randn(Nr,1) + 1j*randn(Nr,1);
noise = noise / sqrt(2*k*snr);
rx_sym = tx_sym+noise;

end