function [pb] = get_ber_bpsk(snr_db)
    snr = 10.^(snr_db/10);
    pb = qfunc(sqrt(2*snr));
end