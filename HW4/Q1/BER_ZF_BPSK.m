function [ber] = BER_ZF_BPSK(Nt,snr_db)
snr = Nt*10.^(snr_db/10);
ber = 0.5*(1-sqrt(snr./(Nt+snr)));

end