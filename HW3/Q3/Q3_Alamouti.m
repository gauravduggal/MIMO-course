close all
clear all
clc


Ntrials = 10000;
scheme = ["BPSK","QPSK","16-QAM"];

snrvec = 20:-0.5:0;
kvec = [1,2,4];
for kidx = 1:length(kvec)
k = kvec(kidx);
ber_imperfectcsi = zeros(size(snrvec));
parfor snr_idx = 1:length(snrvec)
            Nt = 2;  
            Nr = 2;
            snr = snrvec(snr_idx);
            for trials=1:Ntrials
                H = rand(Nr,Nt)+1j*rand(Nr,Nt);
                H = H/sqrt(2);
                
                %data bits, based on the time channel is stationary
                %assuming data rate = 1MSPS and max doppler in channel = 50 Hz
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
                
                for idx = 1:2*k:Npilots-2*k+1
                    tx_bits_slot1 = data_bits(idx:idx+k-1,1);
                    %pick k bits for time slot 1
                    tx_bits_slot2 = data_bits(idx+k:idx+2*k-1,1);
                    %interpret as gray code and convert to a symbol number
                    tx_symbol_num_slot1 = gray2de(tx_bits_slot1,k);
                    tx_symbol_num_slot2 = gray2de(tx_bits_slot2,k);


                    %convert symbol number to complex sample using constellation
                    tx_symbol_complex_slot1 = desym2complexsym(tx_symbol_num_slot1,k);
                    tx_symbol_complex_slot2 = desym2complexsym(tx_symbol_num_slot2,k);
                    s1 = tx_symbol_complex_slot1;
                    s2 = tx_symbol_complex_slot2;
                    
                    %received samples at rx1 at time slot 1 as summation of transmission from tx1, tx2
                    r11 = H(1,1)*s1+H(1,2)*s2;
                    %received samples at rx1 at time slot 2 as summation of transmission from tx1, tx2
                    r21 = -H(1,1)*conj(s2)+H(1,2)*conj(s1);

                    %received samples at rx2 at time slot 1 as summation of transmission from tx1, tx2
                    r12 = H(2,1)*s1+H(2,2)*s2;
                    %received samples at rx2 at time slot 2 as summation of transmission from tx1, tx2
                    r22 = H(2,1)*-conj(s2)+H(2,2)*conj(s1);
                    

                    r11 = add_noise(r11,snr,k);
                    r21 = add_noise(r21,snr,k);
                    r12 = add_noise(r12,snr,k);
                    r22 = add_noise(r22,snr,k);
                    
                    
                    %channel estimate at rx1 from tx1 and tx2
                    H_est(1,1) = H_est(1,1) + 2/(Npilots/k)*(r11*s1'-r21*s2)/2;
                    H_est(1,2) = H_est(1,2) + 2/(Npilots/k)*(r11*s2'+r21*s1)/2;
                    %channel estimate at rx2 from tx1 and tx2
                    H_est(2,1) = H_est(2,1) + 2/(Npilots/k)*(r12*s1'-r22*s2)/2;
                    H_est(2,2) = H_est(2,2) + 2/(Npilots/k)*(r12*s2'+r22*s1)/2;
                        

                end
%                 H_est = H;
%                 %decode the remaining bits using estimated channel
                rx_bits = zeros(Nbits-Npilots,1);
                tx_data_bits = data_bits(Npilots+1:Nbits,1);
                
%                 H_est_eq_rx1 = create_eq_alamouti_rx([H(1,1),H(1,2)]);
%                 pow_rx1 = abs(H(1,1))^2+abs(H(1,2))^2;
%                 H_est_eq_rx2 = create_eq_alamouti_rx([H(2,1),H(2,2)]);
%                 pow_rx2 = abs(H(2,1))^2+abs(H(2,2))^2;

%                 H_eq_rx1 = create_eq_alamouti_rx([H(1,1),H(1,2)]);
%                 H_eq_rx2 = create_eq_alamouti_rx([H(2,1),H(2,2)]);
                
                % decode data bits using estimated channel
                for idx = Npilots+1:2*k:Nbits-2*k+1
%                     idx
                    tx_bits_slot1 = data_bits(idx:idx+k-1,1);
                    %pick k bits for time slot 1
                    tx_bits_slot2 = data_bits(idx+k:idx+2*k-1,1);
                    %interpret as gray code and convert to a symbol number
                    tx_symbol_num_slot1 = gray2de(tx_bits_slot1,k);
                    tx_symbol_num_slot2 = gray2de(tx_bits_slot2,k);


                    %convert symbol number to complex sample using constellation
                    tx_symbol_complex_slot1 = desym2complexsym(tx_symbol_num_slot1,k);
                    tx_symbol_complex_slot2 = desym2complexsym(tx_symbol_num_slot2,k);
                    s1 = tx_symbol_complex_slot1;
                    s2 = tx_symbol_complex_slot2;
                    
                    %received samples at rx1 at time slot 1 as summation of transmission from tx1, tx2
                    r11 = H(1,1)*s1+H(1,2)*s2;
                    %received samples at rx1 at time slot 2 as summation of transmission from tx1, tx2
                    r21 = -H(1,1)*conj(s2)+H(1,2)*conj(s1);

                    %received samples at rx2 at time slot 1 as summation of transmission from tx1, tx2
                    r12 = H(2,1)*s1+H(2,2)*s2;
                    %received samples at rx2 at time slot 2 as summation of transmission from tx1, tx2
                    r22 = -H(2,1)*conj(s2)+H(2,2)*conj(s1);
                    

                    r11 = add_noise(r11,snr,k);
                    r21 = add_noise(r21,snr,k);
                    r12 = add_noise(r12,snr,k);
                    r22 = add_noise(r22,snr,k);
                    
                    

                    % z11 is decision metric for s1 on rx antenna 1
%                     z11 = (H_est(1,1)'*r11+H_est(1,2)*r21')/(abs(H_est(1,1))^2 + abs(H_est(1,2))^2);
%                     % z11 is decision metric for s1 on rx antenna 2
%                     z12 = (H_est(2,1)'*r12+H_est(2,2)*r22')/(abs(H_est(2,1))^2 + abs(H_est(2,2))^2);
%                     z1 = (z11+z12)/2;
                    z11 = (H_est(1,1)'*r11+H_est(1,2)*r21');
                    z12 = (H_est(2,1)'*r12+H_est(2,2)*r22');
                    z1 = (z11+z12)/((abs(H_est(1,1))^2 + abs(H_est(1,2))^2)+(abs(H_est(2,1))^2 + abs(H_est(2,2))^2));



%                     % z21 is decision metric for s2 on rx antenna 1
%                     z21 = (H_est(1,2)'*r11-H_est(1,1)*r21')/(abs(H_est(1,1))^2 + abs(H_est(1,2))^2);
%                     % z22 is decision metric for s2 on rx antenna 2
%                     z22 = (H_est(2,2)'*r12-H_est(2,1)*r22')/(abs(H_est(2,1))^2 + abs(H_est(2,2))^2);
%                     z2 = (z21+z22)/2;

                    z21 = (H_est(1,2)'*r11-H_est(1,1)*r21');
                    % z22 is decision metric for s2 on rx antenna 2
                    z22 = (H_est(2,2)'*r12-H_est(2,1)*r22');
                    z2 = (z21+z22)/((abs(H_est(1,1))^2 + abs(H_est(1,2))^2)+(abs(H_est(2,1))^2 + abs(H_est(2,2))^2));

                    %decode received symbol to closest symbol number in
                    %constellation
                    rx_symbol_num_s1 = decode_complex2de(z1,k);
                    rx_symbol_num_s2 = decode_complex2de(z2,k);
%                     rx_symbol_num_perfectcsi = decode_complex2de(rx_symbol_complex_perfectcsi,k);


                    %convert decoded symbol number to data bits using gray code
                    decoded_bits_s1 = desym2data(rx_symbol_num_s1,k)';
                    decoded_bits_s2 = desym2data(rx_symbol_num_s2,k)';
%                     decoded_bits_perfectcsi = desym2data(rx_symbol_num_perfectcsi ,k)';
                    rx_bits(idx-Npilots:idx-Npilots+k-1,1) = decoded_bits_s1;
                    rx_bits(idx-Npilots+k:idx-Npilots+2*k-1,1) = decoded_bits_s2;
%                     rx_bits_perfectcsi(idx-Npilots:idx-Npilots+k-1,1) = decoded_bits_perfectcsi;
                end

                ber_imperfectcsi(snr_idx) = ber_imperfectcsi(snr_idx) + 1/Ntrials*(1-sum(rx_bits == tx_data_bits)/(Nbits-Npilots));
            end
end
legend_label = "2x2 Alamouti "+scheme(kidx)+" imperfect csi";
semilogy(snrvec,ber_imperfectcsi,'-o','DisplayName',legend_label,'MarkerIndices',1:3:41)
ylim([10^-5,1])
xlim([0,20])
hold on
grid on
end
xlabel("E_b/N_0 dB")
ylabel("P_b")
legend

snrvec_linear = 10.^((snrvec)/10);
u = sqrt(snrvec_linear./(2+snrvec_linear));
Mr = 2;
temp = 0;
for k = 0:2*Mr-1
        temp = temp + nchoosek(2*Mr-1+k,k)*(0.5*(1+u)).^k;
end

Pb = (0.5*(1-u)).^(2*Mr).*temp;
semilogy(snrvec,Pb,'-o','DisplayName',"BPSK-2-branch-MRC-theor.",'MarkerIndices',1:3:41)