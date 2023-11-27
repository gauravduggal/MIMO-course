function [out_samples,n,rx_frame] = add_noise(sig,snr_db)
    snr = 10^(snr_db/10);
    rx_frame = sig;
    n = 1/sqrt(2*snr)*(randn(size(rx_frame))+1j*randn(size(rx_frame)));
    out_samples = rx_frame + n;
end

