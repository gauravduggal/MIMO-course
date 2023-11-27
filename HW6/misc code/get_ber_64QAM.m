function [pb] = get_ber_64QAM(snr_db)
    snr = 10.^(snr_db/10);
    pb = 7/12*qfunc(sqrt(snr/21)).*(1-7/8*qfunc(sqrt(snr/21)));
end