
clear all
clc
close all


scheme = ["BPSK","QPSK","16-QAM"];
kvec = [1,2,4];
figure;
for kidx = 1:length(kvec)
    k = kvec(kidx);
    Ntrials = 10000;
    snrvec = 20:-0.5:0;
    ber_perfectcsi = zeros(size(snrvec));
    ber_theoretical = zeros(size(snrvec));
    parfor snr_idx = 1:length(snrvec)
        snr = snrvec(snr_idx);
        for trials=1:Ntrials
            Nt = 2;
            Nr = 2;

            H = randn(Nr,Nt)+1j*randn(Nr,Nt);
            H = H/sqrt(2);

            Nbits = 10*k;

            %generate N data bits

            data_bits = floor(2*rand(Nbits,1));

            %First Npilots bits are pilots
            Npilots = 0*k;


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
            H_est=H;
            %                 sum(abs(H-H_est).^2)
            %decode the remaining bits using estimated channel
            rx_bits = zeros(Nbits-Npilots,1);
            rx_bits_perfectcsi = zeros(Nbits-Npilots,1);
            tx_data_bits = data_bits(Npilots+1:Nbits,1);
            
            [V,D] = eig(H_est*H_est');
            [max_eig,max_eig_col] = max(diag(D));
            max_eigenvector = V(:,max_eig_col);
            g = transpose(max_eigenvector);
            v = (g*H_est)';
            a = norm(v);
            v = v/a;
            w = g';
            
            
            for idx = Npilots+1:k:Nbits-k+1
                %pick k bits
                tx_bits = data_bits(idx:idx+k-1,1);
                %interpret as gray code and convert to a symbol number
                tx_symbol_num = gray2de(tx_bits,k);

                %convert symbol number to complex sample using constellation
                tx_symbol_complex = desym2complexsym(tx_symbol_num,k);

                %add channel effect
%                 rx_symbol_complex = H*H'*g'*tx_symbol_complex/a;
				
				rx_symbol_complex = H*v*tx_symbol_complex;

                %add thermal noise at rx
                rx_symbol_complex = add_noise(rx_symbol_complex,snr,k);
                
                
                %create decision metric
                z = g*rx_symbol_complex;
%                 z = z / (g*H_est*v);

%                 z = w'*rx_symbol_complex;
%                 
                z= z*a/(g*(H_est*H_est')*g');
                
%               rx_symbol_complex_perfectcsi = H'*rx_symbol_complex/sum(abs(H).^2);

                %decode received symbol to closest symbol number in
                %constellation
                rx_symbol_num = decode_complex2de(z,k);
%                 rx_symbol_num_perfectcsi = decode_complex2de(rx_symbol_complex_perfectcsi,k);


                %convert decoded symbol number to data bits using gray code
                decoded_bits = desym2data(rx_symbol_num,k)';
%                 decoded_bits_perfectcsi = desym2data(rx_symbol_num_perfectcsi ,k)';
                rx_bits(idx-Npilots:idx-Npilots+k-1,1) = decoded_bits;
%                 rx_bits_perfectcsi(idx-Npilots:idx-Npilots+k-1,1) = decoded_bits_perfectcsi;
            end
            ber_perfectcsi(snr_idx) = ber_perfectcsi(snr_idx) + 1/Ntrials*(1-sum(rx_bits == tx_data_bits)/(Nbits-Npilots));
%             %theoretical BER
%             if k == 1
%                 ber_theoretical(1,snr_idx) = ber_theoretical(1,snr_idx) + 1/Ntrials*qfunc(sqrt(2*max_eig/k*10^(snr/10)));
%             elseif k==2
%                 ber_theoretical(1,snr_idx) = ber_theoretical(1,snr_idx) + 1/Ntrials*qfunc(sqrt(2*max_eig/k*10^(snr/10)));
%             else
%                 ber_theoretical(1,snr_idx) = ber_theoretical(1,snr_idx) + 1/Ntrials*1.5/k*erfc(sqrt(max_eig/k*k*10.^(snr/10)/10));
%             end

        end
    end
    legend_label1 = scheme(kidx)+" perfect csi";
    semilogy(snrvec,ber_perfectcsi,'-x','DisplayName',legend_label1,'MarkerIndices',1:3:41);
    grid on
    hold on
            xlabel("E_b/N_0 dB")
        ylabel("P_b")
%     legend_label2 = scheme(kidx)+" theoretical";
%     semilogy(snrvec,ber_theoretical,'-x','DisplayName',legend_label2,'MarkerIndices',1:3:41);
    
    ylim([10^-5,1])

    legend
end

Pb = size(snrvec);
L =2;
for snridx=1:length(snrvec)
    c = 10^(snrvec(snridx)/10);
    Pb(1,snridx) = BER_Maximal_ratio_transmission_2rx_Ltx(L,c); 
end
semilogy(snrvec,Pb,'-x','DisplayName',"BPSK-Theory",'MarkerIndices',1:3:41);
legend
 ylim([10^-5,1])
 
