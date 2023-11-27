function [pb] = get_ber_16QAM(snr_db)
    snr = 10.^(snr_db/10);
    pb = 3/4*qfunc(sqrt(snr/5)).*(1-3/4*qfunc(sqrt(snr/5)));
end