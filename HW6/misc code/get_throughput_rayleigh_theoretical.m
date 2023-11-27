function [th] = get_throughput_rayleigh_theoretical(kvec,snrvec,target_ber,ofdm_symbols,N)
Ntrials = 10000;

% mod = ["BPSK","QPSK","16-QAM","64QAM"];

th = zeros(length(snrvec),length(kvec));
for kidx = 1:length(kvec)
    k = kvec(kidx);
    for snridx = 1:length(snrvec)
        for trial = 1:Ntrials
            snr_db = snrvec(snridx);
            h = randn(1)+1j*randn(1);
            h = h /sqrt(2);
            snr_sc = abs(h)^2*10^(snr_db/10);
            ber = get_ber_awgn2(snr_sc,k);
            if ber > target_ber
                throughput = 0;
            else
                throughput = k*N*ofdm_symbols;
            end
            th(snridx,kidx) = th(snridx,kidx) +  throughput;
        end
        th(snridx,kidx) = th(snridx,kidx)/Ntrials;
    end
end
end
% legend_str = mod(kidx);
% plot(snrvec,th,"-x","DisplayName",legend_str,LineWidth=1.5)
% hold on
% end
% lgd = legend(Location="best");
% set(lgd,'Interpreter','latex');
% set(lgd,'FontSize',12);