function [pb] = get_ber_qpsk(snr_db)
    snr = 10.^(snr_db/10);
    pb = qfunc(sqrt(snr));
end