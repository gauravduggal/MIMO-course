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
tau_rms = 10e-9;
%number of pilots

snrvec = 35:-1:5;
ber =  zeros(size(snrvec));

% sumber of subcarriers
N = 128;
Dn = dftmtx(N)/sqrt(N);
%cyclic prefix > delay spread in samples
cp_n = 0;
%no of ofdm symbols in one packet assuming block fading i.e. prop to
%number of ofdm symbols in a packet based on coherence time
ofdm_symbols = 10;

k_vs_sc_snr = size(N,length(snrvec));
ber_scs_snr = size(N,length(snrvec));
th = size(length(snrvec));
target_ber = 1e-2;

%% Assume first ofdm symbol is used for channel estimation and the channel is
%% feedback to the transmitter where bit loading is performed for the rest of
%% the ofdm symbols within the coherence time

for snridx = 1:length(snrvec)
    snr_db = snrvec(snridx);
    %initialise highest modulation order on all subcarriers
    k_per_sc = 6*ones(N,1);
    %     ber_per_sc = get_ber_64QAM(snr_db)*ones(N,1);
    ber_per_sc = get_ber_awgn(snr_db,6)*ones(N,1);
    enabled_scs = ones(N,1);
    throughput = (1-target_ber)*N*6*ofdm_symbols;
    avg_ber = sum(k_per_sc(enabled_scs).*ber_per_sc(enabled_scs))/((1-target_ber)*sum(k_per_sc(enabled_scs)));
    while(avg_ber>target_ber)
        %find sc with max ber amongst enabled carriers
        [~,idx] = max(ber_per_sc);
        k_per_sc(idx) = get_lower_modulation_scheme(k_per_sc(idx));
%         ber_per_sc(idx) = get_ber_k(k_per_sc(idx),snr_db);
        ber_per_sc(idx) = get_ber_awgn(snr_db,k_per_sc(idx));
        enabled_scs = ~isnan(ber_per_sc);
        avg_ber = sum(enabled_scs)/N *sum(k_per_sc(enabled_scs).*ber_per_sc(enabled_scs))/sum(k_per_sc(enabled_scs));
        %         avg_ber = sum(enabled_scs)/N *sum(k_per_sc(enabled_scs).*ber_per_sc(enabled_scs))/((1-target_ber)*sum(k_per_sc(enabled_scs)));
    end
    ber_scs_snr(1:N,snridx) = ber_per_sc;    
    k_vs_sc_snr(1:N,snridx) = k_per_sc;
end

k_vs_sc_snr
parfor snridx=1:length(snrvec)
    tic
    snr_db = snrvec(snridx);
    ber_snr = 0;
    throughput_snr = 0;
    for trials=1:Ntrials

        k_per_sc = k_vs_sc_snr(:,snridx);
        enabled_scs = (k_vs_sc_snr(:,snridx)>0);
        %         dataEnc = reshape(dataEnc,[k,ofdm_symbols*N]);
        %         sym_dec = transpose(bit2int(temp,k));
        
        Nbits = sum(k_per_sc);
        dataIn = zeros(Nbits*ofdm_symbols,1);
        %modulate symbols across subcarriers,
        %         sym = qammod(sym_dec,2^k,'UnitAveragePower',true);
        %add cyclic prefix

        %assume symbols are in frequency domain and convert to time domain with N time samples
        tx = zeros(N*ofdm_symbols,1);
        tx_ofdm_symbol_i = zeros(N,ofdm_symbols);
        ctr = 1;
        for sidx = 1:ofdm_symbols
            tx_ofdm_symbol_i = zeros(N,1);
            for scidx = 1:N
                if ~enabled_scs(scidx)
                    continue
                end
                Nbits_scidx = k_per_sc(scidx);
                bits = floor(2*rand(Nbits_scidx,1));
                sym_dec = transpose(bit2int(bits,Nbits_scidx));
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
        temp = transpose(temp);

        cp = temp(:,N-cp_n+1:N);
        temp2 = [cp,temp];
        %dont ask why is there a transpose, stupid matlab syntax
        tx_cp = reshape(transpose(temp2),[(N+cp_n)*ofdm_symbols,1]);

        % channel is only awgn
        rx_cp = add_noise_td2(tx_cp,snr_db,N,ofdm_symbols);
        %         rx_cp = tx_cp;

        temp = zeros(N,ofdm_symbols);
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


            %pick received symbols (pilots + data)
            received_data_symbols = temp(1:N,sidx);

            %equalisation with estimated channel

            for scidx = 1:N
                %demodulate using modulation scheme
                if ~enabled_scs(scidx)
                    continue
                end
                Nbits_per_sc = k_per_sc(scidx);
                s_est_dec = qamdemod(received_data_symbols(scidx),2^Nbits_per_sc,'UnitAveragePower',true);
                bits = int2bit(transpose(s_est_dec),Nbits_per_sc);
                dataOut(ctr:ctr+Nbits_per_sc-1,1) = bits;
                ctr = ctr + Nbits_per_sc;
            end
        end
        %             bit_error_current_packet = 1 - 1/(k*N)* sum(sum(int2bit(s_est_dec,k)==int2bit(sym_dec(packet_idx,1),k)));
        %             ber_snr = ber_snr +  1/(packets*Ntrials)*bit_error_current_packet;
        ber_trial = (sum(dataIn~=dataOut)/(Nbits*ofdm_symbols));
        ber_snr = ber_snr + 1/Ntrials * ber_trial;
        throughput_snr = throughput_snr + 1/Ntrials*(1-ber_trial)*Nbits*ofdm_symbols;
    end
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
title("OFDM bit loading -  AWGN")
xlabel("SNR per subcarrier ")
ylabel("BER")



kvec = [1,2,4,6];
mod = ["bpsk","qpsk","16-QAM","64-QAM"];
for kidx = 1:length(kvec)
    %     theory_ber1 = get_ber_k(kvec(kidx),snrvec);
    theory_ber = get_ber_awgn(snrvec,kvec(kidx));
    legend_str = "Theory:mod:"+mod(kidx);
    %     semilogy(snrvec,theory_ber1,"-x","DisplayName",legend_str,LineWidth=1.5);
    hold on
    semilogy(snrvec,theory_ber,"-x","DisplayName",legend_str,LineWidth=1.5);
end
legend_str = "Target BER";
yline(target_ber,"-x","DisplayName",legend_str,LineWidth=1.5)
ylim([1e-4,1])
xlim([5,35])
lgd = legend(Location="best");
set(lgd,'Interpreter','latex');
set(lgd,'FontSize',12);

figure
legend_str = "Bit Loading - AWGN";
plot(snrvec,th,"-x","DisplayName",legend_str,LineWidth=1.5);
hold on
xlabel("snr per subcarrier")
ylabel("throughput")
legend_str = "BPSK";
plot(snrvec,get_throughput_awgn(snrvec,N,ofdm_symbols,1,target_ber),"-x","DisplayName",legend_str,LineWidth=1.5)
legend_str = "QPSK";
plot(snrvec,get_throughput_awgn(snrvec,N,ofdm_symbols,2,target_ber),"-x","DisplayName",legend_str,LineWidth=1.5)
legend_str = "16-QAM";
plot(snrvec,get_throughput_awgn(snrvec,N,ofdm_symbols,4,target_ber),"-x","DisplayName",legend_str,LineWidth=1.5)
legend_str = "64-QAM";
plot(snrvec,get_throughput_awgn(snrvec,N,ofdm_symbols,6,target_ber),"-x","DisplayName",legend_str,LineWidth=1.5)
lgd = legend(Location="best");
set(lgd,'Interpreter','latex');
set(lgd,'FontSize',12);
grid on
title("Throughput")

figure;
imagesc(snrvec,1:N,k_vs_sc_snr)
xlabel("SNR per subcarrier");
ylabel("Subcarrier no.");
title("Bit-loading on subcarriers vs SNR ")
grid on
colorbar('Ticks',[1,2,4,6],...
         'TickLabels',{'BPSK','QPSK','16-QAM','64-QAM'})