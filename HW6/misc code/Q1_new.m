clear all
clc
close all
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%sampling frequency
fs = 100e6;
ts=1/fs;
Ntrials = 100;
%delay spread
tau_rms = 20e-9;

snrvec = 35:-1:5;
ber =  zeros(size(snrvec));

% sumber of subcarriers
N = 128;
Dn = dftmtx(N)/sqrt(N);
%cyclic prefix > delay spread in samples i.e. L
cp_n = 17;
L = 8;
%no of ofdm symbols in one packet assuming block fading i.e. prop to
%number of ofdm symbols in a packet based on coherence time
ofdm_symbols = 10;
k = 2;
npvec = [4,8,16,24,32,48,64];
for npidx = 1:length(npvec)
    Np = npvec(npidx);
    pilot_spacing = floor(N/Np);
    %pilot sub-carrier indices
    pilots_sc_idx = round(linspace(1,N,Np));
    npvec(npidx) = length(pilots_sc_idx);
    Np = npvec(npidx);
    %data sub-carrier indices
    data_sc_idx = 1:N;
    data_sc_idx(pilots_sc_idx) = [];
    Ns = length(data_sc_idx);


    %dft matrices
    QNp = dftmtx(Np)/sqrt(Np);
    QNpL = QNp(:,1:Np);

    parfor snridx=1:length(snrvec)
        tic
        snr_db = snrvec(snridx);
        ber_snr = 0;
        throughput_snr = 0;
        for trials=1:Ntrials
            if mod(trials,100)==0
                trials
            end
            h = transpose(get_freq_selective_channel(tau_rms,fs,N));
            h = h(1:Np);
            h = h./sqrt(sum(abs(h).^2));
            H = fft(h,N)/sqrt(N);
            k_per_sc = k*ones(N,1);
            enabled_scs = (k_per_sc>0);

            Nbits = sum(k_per_sc);
            dataIn = zeros(Nbits*ofdm_symbols,1);

            %assume symbols are in frequency domain and convert to time domain
            %with N time samples / ofdm symbol
            tx = zeros(N*ofdm_symbols,1);
            tx_ofdm_symbol_i = zeros(N,ofdm_symbols);
            tx_symbols = zeros(N,ofdm_symbols);
            ctr = 1;
            pilots = zeros(Np,ofdm_symbols);
            tx_data_symbols = zeros(Ns,ofdm_symbols);
            for sidx = 1:ofdm_symbols
                %ofdm symbol without cp contains N times amples
                pilots_ctr = 1;
                data_symbols_ctr = 1;
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
                    if sum(pilots_sc_idx == scidx)
                        pilots(pilots_ctr,sidx) = sym_scidx;
                        pilots_ctr = pilots_ctr + 1;
                    else
                        tx_data_symbols(data_symbols_ctr,sidx) = sym_dec;
                        data_symbols_ctr = data_symbols_ctr + 1;
                    end
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
%             tx_cp = 1/sqrt((N+cp_n)/(N))*tx_cp;
            %         rx_cp = tx_cp;
            rx_cp = conv(tx_cp,h);

            %        rx_cp = tx_cp;
            rx_cp = add_noise_td2(rx_cp,snr_db,N,ofdm_symbols);

            temp = zeros(N,ofdm_symbols);
            rx_symbols = zeros(N, ofdm_symbols);
            dataOut = zeros(size(dataIn));
            ctr = 1;
            rx_data_symbols = zeros(Ns,ofdm_symbols);
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


                ep = temp(pilots_sc_idx,sidx);
                h_est = (QNpL'*QNpL)^-1*QNpL'*inv(diag(pilots(1:Np,sidx)))*ep/sqrt(Np);
                %freq domain of estimated channel
                H_est = fft(h_est,N)/sqrt(N);

                %equalise channel using estimated channel
                %             received_data_symbols = temp(1:N,sidx);
                received_data_symbols = temp(1:N,sidx)./(sqrt(N)*H_est);

                data_symbols_ctr = 1;
                for scidx = 1:N
                    %demodulate using modulation scheme
                    if ~enabled_scs(scidx)
                        continue
                    end
                    Nbits_per_sc = k_per_sc(scidx);
                    s_est_dec = qamdemod(received_data_symbols(scidx),2^Nbits_per_sc,'UnitAveragePower',true);
                    rx_symbols(scidx,sidx) = s_est_dec;
                    if sum(pilots_sc_idx == scidx)==0
                        rx_data_symbols(data_symbols_ctr,sidx) = s_est_dec;
                        data_symbols_ctr = data_symbols_ctr + 1;
                    end
                    bits = int2bit(transpose(s_est_dec),Nbits_per_sc);
                    dataOut(ctr:ctr+Nbits_per_sc-1,1) = bits;
                    ctr = ctr + Nbits_per_sc;
                end
            end
            %             bit_error_current_packet = 1 - 1/(k*N)* sum(sum(int2bit(s_est_dec,k)==int2bit(sym_dec(packet_idx,1),k)));
            %             ber_snr = ber_snr +  1/(packets*Ntrials)*bit_error_current_packet;
            %ber based on symbol error rate (offers worse performance than
            %theoretical)
%             ber_trial = (sum((sum(tx_data_symbols~=rx_data_symbols)))/(k*Ns*ofdm_symbols));
            
            ber_trial = N/sum(enabled_scs)*(sum(dataIn~=dataOut)/(Nbits*ofdm_symbols));
%             ber_trial = sum(sum((int2bit(tx_data_symbols,2^k) ~= int2bit(rx_data_symbols,2^k))))/(k*Ns*ofdm_symbols);
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

    legend_str ="Npilots-"+string(Np)+" Delay Spread-"+string(tau_rms* 1e9)+"ns Cyclic prefix-"+string(cp_n*ts*1e9)+"ns";
    semilogy(snrvec,ber,"-o","DisplayName",legend_str,LineWidth=1.5);
    % plot(snrvec,mse_channel_estimation,"-x","DisplayName",legend_str,LineWidth=1.5)
    % ylim([1e-5 1])
    % xlim([0,20])
    hold on
    grid on
    title("OFDM Channel Estimation")
    xlabel("SNR per subcarrier")
    ylabel("BER")

end

theory_ber = get_ber_rayleigh(snrvec,k);

legend_str = "QPSK mod Theory";
semilogy(snrvec,theory_ber,"-x","DisplayName",legend_str,LineWidth=1.5);

ylim([1e-4,1])
xlim([5,35])
lgd = legend(Location="best");
set(lgd,'Interpreter','latex');
set(lgd,'FontSize',12);


