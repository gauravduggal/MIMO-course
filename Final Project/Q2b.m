clear all
clc
close all
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');


snrvec = 30:-2:0;
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
Blocks = 500;% size(h,4);
%use bchnumerr(511) and pick a value from second column, last column is
%number of errors
Ndatabits = 250;
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
%         snr_sc_db = 10*log10(N*abs(H).^2*snr);
        %initialise highest modulation order on all subcarriers
%         [k_per_sc,avg_ber,Ne,enabled_scs] = get_bitloading(10.^(snr_sc_db/10),target_ber);
%         if Ne < 1
%             continue
%         end
        k_per_sc = 4*ones(N,1);
        %Number of data bits per ofdm symbol
        Nbits = sum(k_per_sc);
        
        %bchnumerr(n) to get possible total tx bits (n), data bits (k)
        % and error correcting capability
        %n = 511 = s^M-1 i.e. M = 9
        dataIn = floor(2*rand(Ndatabits,1));
        dataInEn = transpose(BCH_encode(dataIn,9,Ndatabits));
        %append dummy bit to make it 512 total bits
        dataInEn = [0; dataInEn];
        temp = bit2int(reshape(dataInEn,[2,256]),2);
        %transmit decimal symbols in slot1 on tx1
        s1 =  qammod(temp(1:128),2^2,'UnitAveragePower',true);
        %transmit decimal symbols in slot2 on tx1
        s2 = qammod(temp(129:256),2^2,'UnitAveragePower',true);
        %-s2* on slot1 tx2, s1* on slot2 tx2
        
        
        tx1 = ifft(s1,N)*sqrt(N);
        tx1_cp = [tx1(:,N-cp_n+1:N), tx1];

        tx2 = ifft(s2,N)*sqrt(N);
        tx2_cp = [tx2(:,N-cp_n+1:N), tx2];

        tx1 = reshape(transpose(tx1_cp),[N+cp_n,1]); 
        tx2 = reshape(transpose(tx2_cp),[N+cp_n,1]);

        %effect of the channel
        rx1 = conv(tx1,h11)+conv(tx2,h12);
        rx2 = conv(tx1,h21)+conv(tx2,h22);
        %receiver noise
        rx1 = add_noise(rx1, snr_db);
        rx2 = add_noise(rx2, snr_db);
        %remove cp and matched filter
        x1 = fft(rx1(cp_n+1:cp_n+1+N),N)/sqrt(N);
        x2 = fft(rx2(cp_n+1:cp_n+1+N),N)/sqrt(N);
        
        H = [diag(H11),diag(H12);diag(H21),diag(H22)];
        r = [x1;x2];
%         temp = (H'*H)^-1*H'*r;
%         I = diag([1./snr_sc_tx1_lin;1./snr_sc_tx2_lin]);
%         temp = (H'*H+I)^-1*H'*r;
        
        temp = sqrt(2)*(H'*H+2/snr*eye(2*N))^-1*H'*r;

        s1_noisy = temp(1:N);
        s2_noisy = temp(N+1:2*N);

        s1_est_dec = qamdemod(s1_noisy,4,'UnitAveragePower',true);
        bits_s1 = int2bit(transpose(s1_est_dec),2);
        bits_s1 = reshape(bits_s1,[256,1]);
        s2_est_dec = qamdemod(s2_noisy,4,'UnitAveragePower',true);
        bits_s2 = int2bit(transpose(s2_est_dec),2);
        bits_s2 = reshape(bits_s2,[256,1]);
        %drop first data bit
        dataOutEn = [bits_s1(2:end);bits_s2];
        dataOut = transpose(BCH_decode(dataOutEn,9,Ndatabits));
        ber_block = sum(dataIn~=dataOut)/Ndatabits;
        if ber_block>0
            throughput = 0;
        else
            throughput = Ndatabits/(TS);
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
legend_str ="MMSE 2x2 MIMO OFDM with BCH coding- Delay Spread:"+string(tau_rms*1e6)+"us";
semilogy(snrvec,ber,"-o","DisplayName",legend_str,LineWidth=1.5);
ylim([1e-5 1])
xlim([0,20])
hold on
grid on
title("PER")
xlabel("snr per subcarrier")
ylabel("PER")



ylim([1e-4,1])
lgd = legend(Location="best");
set(lgd,'Interpreter','latex');
set(lgd,'FontSize',10);


figure
legend_str ="MMSE 2x2 MIMO OFDM with BCH coding- Delay Spread:"+string(tau_rms*1e6)+"us";
plot(snrvec,th,"-x","DisplayName",legend_str,LineWidth=1.5);
hold on
xlabel("snr per subcarrier")
ylabel("throughput")

kvec = 2;
mod = "qpsk";
th = get_throughput_rayleigh_theoretical(kvec,snrvec,bchnumerr(511,Ndatabits)/511,1,N,TS);
legend_str = "Theory:mod:"+mod;
plot(snrvec,th,"-x","DisplayName",legend_str,LineWidth=1.5);


ylim([1 2e6])
lgd = legend(Location="best");
set(lgd,'Interpreter','latex');
set(lgd,'FontSize',10);
grid on
title("Throughput")



