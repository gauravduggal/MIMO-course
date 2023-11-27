
clear all
clc
close all

Lvec = [1,2,4];
scheme = ["BPSK","QPSK","16-QAM"];
kvec = [1,2,4];
for kidx = 1:length(kvec)
    k = kvec(kidx);
    figure;
    for Lidx=1:length(Lvec)
        Nr = Lvec(Lidx);
        
        Ntrials = 4000;
        snrvec = 20:-0.5:0;
        ber_imperfectcsi = zeros(size(snrvec));
        ber_perfectcsi = zeros(size(snrvec));
        parfor snr_idx = 1:length(snrvec)
            tic
            snr = snrvec(snr_idx)
            for trials=1:Ntrials
                H = randn(Nr,1)+1j*randn(Nr,1);
                H = H./sqrt(2);
        
                %data bits, based on the time channel is stationary
                %assuming data rate = 1MSPS and max doppler in channel = 50 Hz
                Nbits = 500*k;

                %generate N data bits

                data_bits = floor(2*rand(Nbits,1));

                %First Npilots bits are pilots
                Npilots = 100*k;

                Nt = 1;


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
            toc
        end
        legend_label1 = string(scheme(kidx))+" L:"+string(Lvec(Lidx))+" K:"+string(kvec(kidx)+" imperfect csi");
        semilogy(snrvec,ber_imperfectcsi,'-x','DisplayName',legend_label1,'MarkerIndices',1:3:41);
        grid on
        hold on
        legend_label2 = string(scheme(kidx))+" L:"+string(Lvec(Lidx))+" K:"+string(kvec(kidx)+" perfect csi");
        semilogy(snrvec,ber_perfectcsi,'-o','DisplayName',legend_label2,'MarkerIndices',1:3:41)
        ylim([10^-5,1])
        xlabel("E_b/N_0 dB")
        ylabel("P_b")
    end 
    legend
end

% legend("L=1, BPSK","L=2, BPSK","L=4, BPSK","L=1, QPSK","L=2, QPSK","L=4, QPSK","L=1, 16QAM","L=2, 16QAM","L=4, 16QAM")