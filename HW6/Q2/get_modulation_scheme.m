function [k] = get_modulation_scheme(received_snr_db,bit_error_rate_target)
    pb_bpsk = get_ber_bpsk(received_snr_db);
    pb_qpsk = get_ber_qpsk(received_snr_db);
    pb_16QAM = get_ber_16QAM(received_snr_db);
    pb_64QAM = get_ber_64QAM(received_snr_db);
    kvec = [0,1,2,4,6];
    idx = find ([pb_bpsk,pb_qpsk,pb_16QAM,pb_64QAM] > bit_error_rate_target,1);
    if ~isempty(idx)
        k = kvec(idx);
    else
        k = 6;
    end
end