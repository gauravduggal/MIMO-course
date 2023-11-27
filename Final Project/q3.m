clear all
clc
close all
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');


snrvec = 35:-2:0;
ber =  zeros(size(snrvec));

% sumber of subcarriers
N = 128;
%cyclic prefix > delay spread in samples i.e. L
cp_n = 40;

k_vs_sc_snr = zeros(N,length(snrvec));

avg_ber_snr = zeros(1,length(snrvec));
enabled_scs_snr = zeros(N,length(snrvec));
th = size(length(snrvec));

load("h.mat");
Blocks = size(h,4);
%use bchnumerr(511) and pick a value from second column, last column is
%number of errors
Ndatabits = 250;
npvec = [8,12,16,20,24,28,32,36,40,50,60];
% npvec=[16,32];
for npidx = 1:length(npvec)
Npilots = npvec(npidx);
pilot_scs_tx1_idx = floor(linspace(2,N-1,Npilots));
pilot_scs_tx2_idx = floor(linspace(2,N-1,Npilots))+1;
data_scs = 1:N;
data_scs(sort([pilot_scs_tx1_idx,pilot_scs_tx2_idx]))=[];

%dft matrices
QNp = dftmtx(Npilots)/sqrt(Npilots);
QNpL = QNp(:,1:Npilots);
Ndata = 2*N-4*Npilots;

parfor snridx=1:length(snrvec)
    tic
    snr_db = snrvec(snridx);
    ber_snr = 0;
    throughput_snr = 0;
    for bidx=1:Blocks

        h_block = h(:,:,:,bidx);
        h11 = squeeze(h_block(1,1,:));
        h12 = squeeze(h_block(1,2,:));
        h21 = squeeze(h_block(2,1,:));
        h22 = squeeze(h_block(2,2,:));

        H11 = fft(h11,N)/sqrt(N);
        H12 = fft(h12,N)/sqrt(N);
        H21 = fft(h21,N)/sqrt(N);
        H22 = fft(h22,N)/sqrt(N);
        %power normalisation to ensure each subcarrier gets unit power
        H11 = sqrt(N)*H11;
        H12 = sqrt(N)*H12;
        H21 = sqrt(N)*H21;
        H22 = sqrt(N)*H22;
        snr = 10^(snr_db/10);

        snr_sc_tx1_lin = abs(H11).^2*snr;
        snr_sc_tx1_db = 10*log10(snr_sc_tx1_lin);
        snr_sc_tx2_lin = abs(H12).^2*snr;
        snr_sc_tx2_db = 10*log10(snr_sc_tx2_lin);
        %initialise highest modulation order on all subcarriers
%         [k_per_sc_tx,avg_ber_tx,Ne_tx,enabled_scs_tx,~] = get_bitloading([10.^(snr_sc_tx1_db/10); 10.^(snr_sc_tx2_db/10)],0.1);
%         k_per_sc_tx1 = k_per_sc_tx(1:N);
%         k_per_sc_tx2 = k_per_sc_tx(N+1:2*N);
        [k_per_sc_tx1,avg_ber_tx1,Ne_tx1,enabled_scs_tx1,~] = get_bitloading(10.^(snr_sc_tx1_db/10),0.1,data_scs);
        [k_per_sc_tx2,avg_ber_tx2,Ne_tx2,enabled_scs_tx2,~] = get_bitloading(10.^(snr_sc_tx2_db/10),0.1,data_scs);

%         k_per_sc_tx1 = 2*ones(N,1);
%         k_per_sc_tx2 = k_per_sc_tx1;
        Ndatabits = sum(k_per_sc_tx1(data_scs))+sum(k_per_sc_tx2(data_scs));


        %packet contains two ofdm symbols each with N subcarriers. This is
        %per tx antenna so one packet contains N*4 symbols. The first ofdm
        %symbol contains the pilots for channel estimation
        [tx1_complex_sym,tx2_complex_sym,tx_data_bits,tx_data_symbols] = tx_place_data_sc(k_per_sc_tx1,k_per_sc_tx2,data_scs);

        %place pilots - assumed to be QPSK for simplicity
        [tx1_complex_sym,tx2_complex_sym, pilots_tx1, pilots_tx2] = tx_place_pilots(tx1_complex_sym,tx2_complex_sym, pilot_scs_tx1_idx, pilot_scs_tx2_idx);

        

        tx1 = ifft(tx1_complex_sym,N)*sqrt(N);
        tx1_cp = [tx1(N-cp_n+1:N,:); tx1];

        tx2 = ifft(tx2_complex_sym,N)*sqrt(N);
        tx2_cp = [tx2(N-cp_n+1:N,:); tx2];

        tx1 = reshape(transpose(tx1_cp),[N+cp_n,1]);
        tx2 = reshape(transpose(tx2_cp),[N+cp_n,1]);

        %effect of the channel
%         rx1 = tx1;
%         rx2 = tx2;
        rx1 = conv(tx1,h11)+conv(tx2,h12);
        rx2 = conv(tx1,h21)+conv(tx2,h22);
        %receiver noise
        rx1 = add_noise(rx1, snr_db);
        rx2 = add_noise(rx2, snr_db);
        %remove cp and matched filter
        x1 = fft(rx1(cp_n+1:cp_n+N),N)/sqrt(N);
        x2 = fft(rx2(cp_n+1:cp_n+N),N)/sqrt(N);


        ep_rx1_tx1_noisy_pilots = x1(pilot_scs_tx1_idx);
        ep_rx1_tx2_noisy_pilots = x1(pilot_scs_tx2_idx);
        ep_rx2_tx1_noisy_pilots = x2(pilot_scs_tx1_idx);
        ep_rx2_tx2_noisy_pilots = x2(pilot_scs_tx2_idx);

        h11_est = (QNpL'*QNpL)^-1*QNpL'*inv(diag(pilots_tx1))*ep_rx1_tx1_noisy_pilots/sqrt(Npilots);
        h12_est = (QNpL'*QNpL)^-1*QNpL'*inv(diag(pilots_tx2))*ep_rx1_tx2_noisy_pilots/sqrt(Npilots);
        h21_est = (QNpL'*QNpL)^-1*QNpL'*inv(diag(pilots_tx1))*ep_rx2_tx1_noisy_pilots/sqrt(Npilots);
        h22_est = (QNpL'*QNpL)^-1*QNpL'*inv(diag(pilots_tx2))*ep_rx2_tx2_noisy_pilots/sqrt(Npilots);
        
        H11_est = fft(h11_est,N)/sqrt(N);
        H12_est = fft(h12_est,N)/sqrt(N);
        H21_est = fft(h21_est,N)/sqrt(N);
        H22_est = fft(h22_est,N)/sqrt(N);
        %power normalisation to ensure each subcarrier gets unit power
        H11_est = sqrt(N)*H11_est;
        H12_est = sqrt(N)*H12_est;
        H21_est = sqrt(N)*H21_est;
        H22_est = sqrt(N)*H22_est;

        H_est = [diag(H11_est),diag(H12_est);diag(H21_est),diag(H22_est)];

        H = [diag(H11),diag(H12);diag(H21),diag(H22)];
        r = [x1;x2];
%         temp = r;
%         temp = (H'*H)^-1*H'*r;
%         temp = (H_est'*H_est)^-1*H_est'*r;
        I = diag([1./snr_sc_tx1_lin;1./snr_sc_tx2_lin]);
       temp = sqrt(2)*(H_est'*H_est+2/snr*eye(2*N))^-1*H_est'*r;
%         temp = (H'*H+1/snr*eye(2*N))^-1*H'*r;
%       
        rx1_equalized_noisy = temp(1:N);
        rx2_equalized_noisy = temp(N+1:2*N);
        

        [rx_data_bits,rx_data_symbols] = rx_get_data_sc(rx1_equalized_noisy,rx2_equalized_noisy, k_per_sc_tx1,k_per_sc_tx2,data_scs);
        
        ber_block = sum(tx_data_bits~=rx_data_bits)/length(tx_data_bits);
        if ber_block>0.25
            throughput = 0;
        else
            throughput = length(tx_data_bits)/(TS);
        end
        throughput_snr = throughput_snr + 1/Blocks *throughput;
        ber_snr = ber_snr + 1/Blocks* ber_block;

    end
    %     deep_fade_instantiations
    ber(snridx) = ber_snr;
    %     th(snridx) = throughput_snr/Ntrials;
    th(snridx) = throughput_snr;
    toc
end
% legend_str ="Alamouti 2x2 MIMO OFDM: Delay Spread:"+string(tau_rms*1e9)+"ns"+" Pilots"+string(Npilots);
% semilogy(snrvec,ber,"-o","DisplayName",legend_str,LineWidth=1.5);
% ylim([1e-5 1])
% xlim([0,20])
% hold on
% grid on
% title("OFDM ")
% xlabel("snr per subcarrier")
% ylabel("PER")
% 
% 
% 
% ylim([1e-4,1])
% lgd = legend(Location="best");
% set(lgd,'Interpreter','latex');
% set(lgd,'FontSize',10);
% 


legend_str = " Np:"+string(Npilots);
plot(snrvec,th,"-x","DisplayName",legend_str,LineWidth=1.5);
hold on
xlabel("snr per subcarrier")
ylabel("throughput")
end
% kvec = 2;
% mod = "qpsk";
% th = get_throughput_rayleigh_theoretical(kvec,snrvec,bchnumerr(511,Ndatabits)/511,2,N,2*TS);
% legend_str = "Theory:mod:"+mod;
% semilogy(snrvec,th,"-x","DisplayName",legend_str,LineWidth=1.5);


ylim([0 12e6])
lgd = legend(Location="best");
set(lgd,'Interpreter','latex');
set(lgd,'FontSize',10);
grid on
title("Throughput-2X2 MIMO-OFDM spatial multiplexing + BitLoading")



