
clear all
clc
close all


scheme = ["BPSK","QPSK","16-QAM"];
kvec = [1,2,4];
figure;
for kidx = 1:length(kvec)
    k = kvec(kidx);
    Ntrials = 20000;
    snrvec = 20:-0.5:0;
    ber_imperfectcsi = zeros(size(snrvec));
    ber_perfectcsi = zeros(size(snrvec));
    parfor snr_idx = 1:length(snrvec)
        snr = snrvec(snr_idx);
        for trials=1:Ntrials
            Nt = 2;
            Nr = 2;

            H = randn(Nr,Nt)+1j*randn(Nr,Nt);
            H = H./sqrt(2);

            Nbits = 500*k;

            %generate N data bits

            data_bits = floor(2*rand(Nbits,1));

            %First Npilots bits are pilots
            Npilots = 100*k;


            %generate channel for 1 packet containing Nbits total bits with the first
            %Npilots pilot bits
            %                 H = randn(Nr,Nt)+1j*randn(Nr,Nt);

            %bits / symbol
            %             k = 1;

            H_est = zeros(Nr,Nt);
            % channel estimation using pilots
            for idx = 1:k:Npilots-k+1
                %pick k bits
                tx_bits = data_bits(idx:idx+k-1,1);
                %interpret as gray code and convert to a symbol number
                tx_symbol_num = gray2de(tx_bits,k);

                %convert symbol number to complex sample using constellation
                tx_symbol_complex = desym2complexsym(tx_symbol_num,k);

                %add channel effect
                rx_symbol_complex = H*tx_symbol_complex;

                %add thermal noise at rx
                rx_symbol_complex = add_noise(rx_symbol_complex,snr,k);

                H_est = H_est + 1/(Npilots/k) * rx_symbol_complex*tx_symbol_complex';

            end
            %                 sum(abs(H-H_est).^2)
            %decode the remaining bits using estimated channel
            rx_bits = zeros(Nbits-Npilots,1);
            rx_bits_perfectcsi = zeros(Nbits-Npilots,1);
            tx_data_bits = data_bits(Npilots+1:Nbits,1);
            %Channel matrix for tx1
            h1 = H_est(:,1);
            %Channel matrix for tx2
            h2 = H_est(:,2);
            
%             [~,D]= eig(h1*h1');
%             eig_h = max(diag(D));
%             C1 = log2(1+(10^(k*snr/10)*eig_h));
%             
%             [~,D]= eig(h2*h2');
%             eig_h = max(diag(D));
%             C2 = log2(1+(10^(k*snr/10)*eig_h));
%             
            %Capacity if using tx1


%             if C2>C1
%                 H_est = h2;
%                 H = H(:,2);
%                 
%             else
%                 H_est = h1;
%                 H = H(:,1);
%             end
%             
            if sum(abs(h2).^2) > sum(abs(h1).^2)
                H_est = h2;
                H=H(:,2);
            else
                H_est = h1;
                H=H(:,1);
            end
            % decode data bits using estimated channel

            for idx = Npilots+1:k:Nbits-k+1
                %pick k bits
                tx_bits = data_bits(idx:idx+k-1,1);
                %interpret as gray code and convert to a symbol number
                tx_symbol_num = gray2de(tx_bits,k);

                %convert symbol number to complex sample using constellation
                tx_symbol_complex = desym2complexsym(tx_symbol_num,k);

                %add channel effect
                rx_symbol_complex = H*tx_symbol_complex;

                %add thermal noise at rx
                rx_symbol_complex = add_noise(rx_symbol_complex,snr,k);

                %remove channel effect using maximal ratio combining using estimated channel
                rx_symbol_complex_est = H_est'*rx_symbol_complex/sum(abs(H_est).^2);

                rx_symbol_complex_perfectcsi = H'*rx_symbol_complex/sum(abs(H).^2);

                %decode received symbol to closest symbol number in
                %constellation
                rx_symbol_num = decode_complex2de(rx_symbol_complex_est,k);
                rx_symbol_num_perfectcsi = decode_complex2de(rx_symbol_complex_perfectcsi,k);


                %convert decoded symbol number to data bits using gray code
                decoded_bits = desym2data(rx_symbol_num,k)';
                decoded_bits_perfectcsi = desym2data(rx_symbol_num_perfectcsi ,k)';
                rx_bits(idx-Npilots:idx-Npilots+k-1,1) = decoded_bits;
                rx_bits_perfectcsi(idx-Npilots:idx-Npilots+k-1,1) = decoded_bits_perfectcsi;
            end
            ber_imperfectcsi(snr_idx) = ber_imperfectcsi(snr_idx) + 1/Ntrials*(1-sum(rx_bits == tx_data_bits)/(Nbits-Npilots));
            ber_perfectcsi(snr_idx) = ber_perfectcsi(snr_idx) + 1/Ntrials*(1-sum(rx_bits_perfectcsi == tx_data_bits)/(Nbits-Npilots));
        end
    end
    legend_label1 = scheme(kidx)+" imperfect csi";
    semilogy(snrvec,ber_imperfectcsi,'-x','DisplayName',legend_label1,'MarkerIndices',1:3:41);
    grid on
    hold on
    legend_label2 = scheme(kidx)+" perfect csi";
%     semilogy(snrvec,ber_perfectcsi,'-o','DisplayName',legend_label2,'MarkerIndices',1:3:41)
    ylim([10^-5,1])

    legend
end
    k = 1;
    Pb = zeros(size(snrvec));
    for snridx = 1:length(snrvec)
        snr = snrvec(snridx);
        gamma = 10^(snr/10);
        temp = 0;
        Nt = 2;
        for k = 1:Nt
            temp = temp + nchoosek(Nt,k)*(-1)^(k-1)*sqrt(gamma/(gamma+k));
        end
        Pb(snridx) = 0.5*(1-temp);
    end
    hold on
%     semilogy(snrvec,Pb,'-o','DisplayName',"Theoretical",'MarkerIndices',1:3:41)

        xlabel("E_b/N_0 dB")
        ylabel("P_b")


for kidx = 1:length(kvec)
    k = kvec(kidx);
    Ntrials = 20000;
    snrvec = 20:-0.5:0;
    ber_imperfectcsi = zeros(size(snrvec));
    ber_perfectcsi = zeros(size(snrvec));
    parfor snr_idx = 1:length(snrvec)
        snr = snrvec(snr_idx);
        for trials=1:Ntrials
            Nt = 2;
            Nr = 2;

            H = randn(Nr,Nt)+1j*randn(Nr,Nt);
            H = H./sqrt(2);

            Nbits = 500*k;

            %generate N data bits

            data_bits = floor(2*rand(Nbits,1));

            %First Npilots bits are pilots
            Npilots = 100*k;


            %generate channel for 1 packet containing Nbits total bits with the first
            %Npilots pilot bits
            %                 H = randn(Nr,Nt)+1j*randn(Nr,Nt);

            %bits / symbol
            %             k = 1;

            H_est = zeros(Nr,Nt);
            % channel estimation using pilots
            for idx = 1:k:Npilots-k+1
                %pick k bits
                tx_bits = data_bits(idx:idx+k-1,1);
                %interpret as gray code and convert to a symbol number
                tx_symbol_num = gray2de(tx_bits,k);

                %convert symbol number to complex sample using constellation
                tx_symbol_complex = desym2complexsym(tx_symbol_num,k);

                %add channel effect
                rx_symbol_complex = H*tx_symbol_complex;

                %add thermal noise at rx
                rx_symbol_complex = add_noise(rx_symbol_complex,snr,k);

                H_est = H_est + 1/(Npilots/k) * rx_symbol_complex*tx_symbol_complex';

            end
            %                 sum(abs(H-H_est).^2)
            %decode the remaining bits using estimated channel
            rx_bits = zeros(Nbits-Npilots,1);
            rx_bits_perfectcsi = zeros(Nbits-Npilots,1);
            tx_data_bits = data_bits(Npilots+1:Nbits,1);
            %Channel matrix for tx1
            h1 = H_est(:,1);
            %Channel matrix for tx2
            h2 = H_est(:,2);
            
%             [~,D]= eig(h1*h1');
%             eig_h = max(diag(D));
%             C1 = log2(1+(10^(k*snr/10)*eig_h));
%             
%             [~,D]= eig(h2*h2');
%             eig_h = max(diag(D));
%             C2 = log2(1+(10^(k*snr/10)*eig_h));
%             
%             %Capacity if using tx1
% %             C1 = log2(1+(10^(k*snr/10)*abs(h1(1,1))^2)/2) + log2(1+(10^(k*snr/10)*abs(h1(2,1))^2)/2);
%             %Capacity if using tx2
% %             C2 = log2(1+(10^(k*snr/10)*abs(h2(1,1))^2)/2) + log2(1+(10^(k*snr/10)*abs(h2(2,1))^2)/2);
%             if C2>C1 && rand(1) > 0.1
%                 H_est = h2;
%                 H = H(:,2);
%                 
%             else
%                 H_est = h1;
%                 H = H(:,1);
%             end
%             
            if sum(abs(h2).^2) > sum(abs(h1).^2) && rand(1)>0.1
                H_est = h2;
                H=H(:,2);
            else
                H_est = h1;
                H=H(:,1);
            end
            % decode data bits using estimated channel

            for idx = Npilots+1:k:Nbits-k+1
                %pick k bits
                tx_bits = data_bits(idx:idx+k-1,1);
                %interpret as gray code and convert to a symbol number
                tx_symbol_num = gray2de(tx_bits,k);

                %convert symbol number to complex sample using constellation
                tx_symbol_complex = desym2complexsym(tx_symbol_num,k);

                %add channel effect
                rx_symbol_complex = H*tx_symbol_complex;

                %add thermal noise at rx
                rx_symbol_complex = add_noise(rx_symbol_complex,snr,k);

                %remove channel effect using maximal ratio combining using estimated channel
                rx_symbol_complex_est = H_est'*rx_symbol_complex/sum(abs(H_est).^2);

                rx_symbol_complex_perfectcsi = H'*rx_symbol_complex/sum(abs(H).^2);

                %decode received symbol to closest symbol number in
                %constellation
                rx_symbol_num = decode_complex2de(rx_symbol_complex_est,k);
                rx_symbol_num_perfectcsi = decode_complex2de(rx_symbol_complex_perfectcsi,k);


                %convert decoded symbol number to data bits using gray code
                decoded_bits = desym2data(rx_symbol_num,k)';
                decoded_bits_perfectcsi = desym2data(rx_symbol_num_perfectcsi ,k)';
                rx_bits(idx-Npilots:idx-Npilots+k-1,1) = decoded_bits;
                rx_bits_perfectcsi(idx-Npilots:idx-Npilots+k-1,1) = decoded_bits_perfectcsi;
            end
            ber_imperfectcsi(snr_idx) = ber_imperfectcsi(snr_idx) + 1/Ntrials*(1-sum(rx_bits == tx_data_bits)/(Nbits-Npilots));
            ber_perfectcsi(snr_idx) = ber_perfectcsi(snr_idx) + 1/Ntrials*(1-sum(rx_bits_perfectcsi == tx_data_bits)/(Nbits-Npilots));
        end
    end
    legend_label1 = scheme(kidx)+"10% feedback error"+" imperfect csi";
    semilogy(snrvec,ber_imperfectcsi,'-x','DisplayName',legend_label1,'MarkerIndices',1:3:41);
    grid on
    hold on
%     legend_label2 = scheme(kidx)+" perfect csi";
%     semilogy(snrvec,ber_perfectcsi,'-o','DisplayName',legend_label2,'MarkerIndices',1:3:41)
    ylim([10^-5,1])

    legend
            xlabel("E_b/N_0 dB")
        ylabel("P_b")
end


Pb = factorial(4*2-1)./(2^(5*2)*factorial(2*2-1)).*(1./10.^(snrvec/10)).^(4);
semilogy(snrvec,Pb,'-o','DisplayName',"BPSK-Theoretical-approx",'MarkerIndices',1:3:41);

