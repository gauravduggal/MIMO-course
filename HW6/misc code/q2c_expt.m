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
tau_rms = 100e-9;
%channel taps
L = 10;
snrvec = 45:-1:5;
ber =  zeros(size(snrvec));

% sumber of subcarriers
N = 128;
Dn = dftmtx(N)/sqrt(N);
%cyclic prefix > delay spread in samples i.e. L
cp_n = 11;
%no of ofdm symbols in one packet assuming block fading i.e. prop to
%number of ofdm symbols in a packet based on coherence time
ofdm_symbols = 10;

k_vs_sc_snr = zeros(N,length(snrvec));
target_ber = 1e-2;
avg_ber_snr = zeros(1,length(snrvec));
enabled_scs_snr = zeros(N,length(snrvec));
th = size(length(snrvec));


for snridx=1:length(snrvec)
    tic
    snr_db = snrvec(snridx);
    ber_snr = 0;
    deep_fade_instantiations = 0;
    throughput_snr = 0;
    for trials=1:Ntrials
%         h = [1, 0.5, 0.35, 0.15, 0.05];
        h = get_freq_selective_channel(tau_rms,fs,L);
%         h = h/sqrt(sum(abs(h).^2));
        H = fft(h,N)/sqrt(N);
%         H = sqrt(N)*H;
        snr = 10^(snr_db/10);
        snr_sc_db = 10*log10(N*abs(H).^2*snr);
        %initialise highest modulation order on all subcarriers
        k_per_sc = 6*ones(N,1);
        theory_ber = transpose(get_ber_awgn(snr_sc_db,6));
        ber_per_sc = theory_ber;
        enabled_scs = ones(N,1);
        %         throughput = target_ber*N*6*ofdm_symbols;
        avg_ber = sum(k_per_sc(enabled_scs).*ber_per_sc(enabled_scs))/sum(k_per_sc(enabled_scs))*sum(enabled_scs)/N;
        while(avg_ber>target_ber)
            %find sc with max ber amongst enabled carriers
            [~,idx] = max(ber_per_sc);
            k_per_sc(idx) = get_lower_modulation_scheme(k_per_sc(idx));
            %             ber_per_sc(idx) = get_ber_k(k_per_sc(idx),snr_sc_db(idx));
            ber_per_sc(idx) = get_ber_awgn(snr_sc_db(idx),k_per_sc(idx));
%             ber_per_sc(idx) = get_ber_rayleigh(snr_db,k_per_sc(idx));
        
            enabled_scs = (ber_per_sc>0);
            avg_ber = sum(enabled_scs)/N* sum(k_per_sc.*ber_per_sc)/sum(k_per_sc);
        end
%         kvec = [0,1,0,0,0];
%         k_per_sc = kvec(round((3*rand(128,1))+1)); 
%         enabled_scs = (k_per_sc>0);
        
        Nbits = sum(k_per_sc);
        dataIn = zeros(Nbits*ofdm_symbols,1);
        %modulate symbols across subcarriers,
        %         sym = qammod(sym_dec,2^k,'UnitAveragePower',true);
        %add cyclic prefix

        %assume symbols are in frequency domain and convert to time domain
        %with N time samples / ofdm symbol
        tx = zeros(N*ofdm_symbols,1);
        tx_ofdm_symbol_i = zeros(N,ofdm_symbols);
        tx_symbols = zeros(N,ofdm_symbols);
        ctr = 1;
        for sidx = 1:ofdm_symbols
            %ofdm symbol without cp contains N times amples
            for scidx = 1:N
                if ~enabled_scs(scidx)
                    continue
                end
                %pick Nbits_scidx on subcarrier scidx
                Nbits_scidx = k_per_sc(scidx);
                %generate required data bits
                bits = floor(2*rand(Nbits_scidx,1));
                %modulate using scheme according to SNR
                sym_dec = bit2int(bits,Nbits_scidx);
                tx_symbols(scidx,sidx) = sym_dec;
                sym_scidx = qammod(sym_dec,2^Nbits_scidx,'UnitAveragePower',true);
                tx_ofdm_symbol_i(scidx,sidx) = sym_scidx;
                dataIn(ctr:ctr+Nbits_scidx-1) = bits;
                ctr = ctr + Nbits_scidx;
            end

            st_idx = (sidx-1)*N+1;
            en_idx = st_idx+N-1;
            %modulate tx symbols
            tx(st_idx:en_idx,1) = ifft(tx_ofdm_symbol_i(:,sidx),N)*sqrt(N);
        end
        %yucky matlab syntax section to add cyclic prefix
        temp = reshape(tx,[N,ofdm_symbols]);
%         temp = transpose(temp);

        cp = temp(N-cp_n+1:N,:);
        temp2 = [cp;temp];
        tx_cp = reshape(temp2,[(N+cp_n)*ofdm_symbols,1]);
        %compensate for power lost in addition of cp
%         tx_cp = 1/sqrt((N+cp_n)/(N))*tx_cp;
        
        rx_cp = conv(tx_cp,h);

        %        rx_cp = tx_cp;
        rx_cp = add_noise_td2(rx_cp,snr_db,N,ofdm_symbols);
        
        temp = zeros(N,ofdm_symbols);
        rx_symbols = zeros(N, ofdm_symbols);
        dataOut = zeros(size(dataIn));
        ctr = 1;
        for sidx = 1:ofdm_symbols
            %strip cyclic prefix from every OFDM symbol and the tail from the
            %last OFDM symbol to get integration length for each OFDM symbol
            st_idx = (sidx-1)*N+(sidx)*cp_n+1;
            en_idx = st_idx+N-1;
            %get integration length
            temp(1:N,sidx) = rx_cp(st_idx:en_idx,1);
            %matched filter  to convert to received distorted symbols in freq domain
            temp(1:N,sidx) = fft(temp(1:N,sidx),N)/(sqrt(N));
            %                 temp(sidx,1:N) = Dn*transpose(temp(sidx,1:N));

            %equalise channel assuming perfect channel estimation
%             received_data_symbols = temp(1:N,sidx);
            received_data_symbols = temp(1:N,sidx)./(sqrt(N)*transpose(H));
            

            for scidx = 1:N
                %demodulate using modulation scheme
                if ~enabled_scs(scidx)
                    continue
                end
                Nbits_per_sc = k_per_sc(scidx);
                s_est_dec = qamdemod(received_data_symbols(scidx),2^Nbits_per_sc,'UnitAveragePower',true);
                rx_symbols(scidx,sidx) = s_est_dec;
                bits = int2bit(transpose(s_est_dec),Nbits_per_sc);
                dataOut(ctr:ctr+Nbits_per_sc-1,1) = bits;
                ctr = ctr + Nbits_per_sc;
            end
        end
        %             bit_error_current_packet = 1 - 1/(k*N)* sum(sum(int2bit(s_est_dec,k)==int2bit(sym_dec(packet_idx,1),k)));
        %             ber_snr = ber_snr +  1/(packets*Ntrials)*bit_error_current_packet;
        ber_trial = N/sum(enabled_scs)*(sum(dataIn~=dataOut)/(Nbits*ofdm_symbols));
       
%         ber_trial = avg_ber;
        ber_snr = ber_snr + 1/Ntrials *  ber_trial;
%         throughput_snr = throughput_snr + 1/Ntrials*(1-ber_trial)*sum(k_per_sc(enabled_scs))*sum(enabled_scs)/N*ofdm_symbols; 
        throughput_snr = throughput_snr + 1/Ntrials*(1-ber_trial)*(Nbits*ofdm_symbols);

    end
%     deep_fade_instantiations
    ber(snridx) = ber_snr;
    th(snridx) = throughput_snr; 
    toc
end
legend_str ="OFDM Bit-loading";
semilogy(snrvec,ber,"-o","DisplayName",legend_str,LineWidth=1.5);
% plot(snrvec,mse_channel_estimation,"-x","DisplayName",legend_str,LineWidth=1.5)
% ylim([1e-5 1])
% xlim([0,20])
hold on
grid on
title("OFDM bit loading -  Freq Selective")
xlabel("SNR (dB)")
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
xlim([5,45])
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

