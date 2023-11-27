function [pb] = get_ber_k(k,snr_db)
switch k
    case 0
        pb = nan;
    case 1
        pb = get_ber_bpsk(snr_db);
    case 2
        pb = get_ber_qpsk(snr_db);
    case 4
        pb = get_ber_16QAM(snr_db);
    case 6
        pb = get_ber_64QAM(snr_db);
    otherwise
        pb = nan;
end
end