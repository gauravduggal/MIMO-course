clear all
clc
close all
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%sampling frequency
fs = 100e6;
ts=1/fs;
Ntrials = 10;
%delay spread
tau_rms = 50e-9;
L=10;
snrvec = 15:-1:10;
ber =  zeros(size(snrvec));

% sumber of subcarriers
N = 128;
Dn = dftmtx(N)/sqrt(N);
%cyclic prefix > delay spread in samples
cp_n = 20;
%no of ofdm symbols in one packet assuming block fading i.e. prop to
%number of ofdm symbols in a packet based on coherence time
ofdm_symbols = 10;

k_vs_sc_snr = zeros(N,length(snrvec));
target_ber = 1e-2;
avg_ber_snr = zeros(1,length(snrvec));
enabled_scs_snr = zeros(N,length(snrvec));
th = zeros(1,length(snrvec));
for snridx = 1:length(snrvec)
        tic
    for trial = 1:Ntrials
        
        snr_db = snrvec(snridx);
        h = get_freq_selective_channel(tau_rms,fs,L);
        H = transpose(fft(h,N)/sqrt(N));
        snr = 10^(snr_db/10);
        %received power per subcarrier assuming total transmit power is N
        snr_sc_db = 10*log10(N*abs(H).^2*snr);
        %initialise highest modulation order on all subcarriers
        [k_per_sc,avg_ber,Ne,enabled_scs,snr_w_null_db] = get_bitloading(10.^(snr_sc_db/10),target_ber);
        %     k_per_sc = 6*ones(N,1);
        %     theory_ber = transpose(get_ber_awgn(snr_sc_db,6));
        %     ber_per_sc = theory_ber;
        %     enabled_scs = ones(N,1);
        %     %         throughput = target_ber*N*6*ofdm_symbols;
        %     avg_ber = sum(enabled_scs)/N * sum(k_per_sc(enabled_scs).*ber_per_sc(enabled_scs))/(sum(k_per_sc(enabled_scs)));
        %     while(avg_ber>target_ber)
        %         %find sc with min ber amongst enabled carriers
        %         [~,idx] = max(ber_per_sc);
        %         k_per_sc(idx) = get_lower_modulation_scheme(k_per_sc(idx));
        %         %             ber_per_sc(idx) = get_ber_k(k_per_sc(idx),snr_sc_db(idx));
        %         ber_per_sc(idx) = get_ber_awgn(snr_sc_db(idx),k_per_sc(idx));
        %         enabled_scs = (ber_per_sc>0);
        %         %             avg_ber = sum(k_per_sc(enabled_scs).*ber_per_sc(enabled_scs))/(mean(k_per_sc)*N);
        %
        %         avg_ber = sum(enabled_scs)/N * sum(k_per_sc(enabled_scs).*ber_per_sc(enabled_scs))/(sum(k_per_sc(enabled_scs)));
        %     end
        k_vs_sc_snr(1:N,snridx) = k_vs_sc_snr(1:N,snridx) + 1/Ntrials*k_per_sc;
        avg_ber_snr(1,snridx) = avg_ber_snr(1,snridx) + 1/Ntrials* avg_ber;
        enabled_scs_snr(:,snridx) = enabled_scs;
        throughput = (1-avg_ber)*ofdm_symbols*sum(k_per_sc);
        th(snridx) = th(snridx) + 1/Ntrials*throughput; 
        
    end
    toc
end








figure;
plot(1:N,snr_sc_db,LineWidth=1.5)
grid on
hold on
xlabel("Subcarrier")
ylabel("SNR over sub-carriers")
scs_idxs = 1:10:N;
text(scs_idxs,snr_sc_db(scs_idxs),string(k_per_sc(scs_idxs))+"bps")
title("Modulation scheme vs Sub-carrier")
% xline(1:N)
figure;
legend_str ="OFDM Bit-loading";
semilogy(snrvec,avg_ber_snr,"-o","DisplayName",legend_str,LineWidth=1.5);
hold on
grid on
title("OFDM bit loading -  Freq Selective")
xlabel("SNR per sub-carrier")
ylabel("BER")



kvec = [1,2,4,6];
mod = ["bpsk","qpsk","16-QAM","64-QAM"];
for kidx = 1:length(kvec)
    theory_ber = get_ber_rayleigh(snrvec,kvec(kidx));
    %     if kidx>1
    %         theory_ber = berfading(10*log10(1/kvec(kidx)*10.^(snrvec/10)),'qam',2^kvec(kidx),1);
    %     else
    %         theory_ber = berfading(10*log10(1/kvec(kidx)*10.^(snrvec/10)),'psk',2^kvec(kidx),1);
    %     end
    %     theory_ber = get_ber_k(kvec(kidx),10*log10(10.^(snrvec/10)));
    legend_str = "Theory:mod:"+mod(kidx);
    semilogy(snrvec,theory_ber,"-x","DisplayName",legend_str,LineWidth=1.5);
    hold on
end
legend_str = "Target BER";
yline(target_ber,"-x","DisplayName",legend_str,LineWidth=1.5)
ylim([1e-4,1])
xlim([10,35])
lgd = legend(Location="best");
set(lgd,'Interpreter','latex');
set(lgd,'FontSize',12);

figure
legend_str = "Bit Loading - Rayleigh Freq. Sel. Fading";
plot(snrvec,th,"-x","DisplayName",legend_str,LineWidth=1.5);
hold on
xlabel("snr")
ylabel("throughput")
legend_str = "BPSK";
plot(snrvec,get_throughput_rayleigh(snrvec,N,ofdm_symbols,1,target_ber),"-x","DisplayName",legend_str,LineWidth=1.5)
legend_str = "QPSK";
plot(snrvec,get_throughput_rayleigh(snrvec,N,ofdm_symbols,2,target_ber),"-x","DisplayName",legend_str,LineWidth=1.5)
legend_str = "16-QAM";
plot(snrvec,get_throughput_rayleigh(snrvec,N,ofdm_symbols,4,target_ber),"-x","DisplayName",legend_str,LineWidth=1.5)
legend_str = "64-QAM";
plot(snrvec,get_throughput_rayleigh(snrvec,N,ofdm_symbols,6,target_ber),"-x","DisplayName",legend_str,LineWidth=1.5)
lgd = legend(Location="best");
set(lgd,'Interpreter','latex');
set(lgd,'FontSize',12);
grid on
title("Throughput")

%  plot(snrvec, avg_ber_snr)
% semilogy(snrvec, avg_ber_snr)
% % ylim([1e-4,1])
