close all
clear all
clc





%number of bits in a packet
N = 100000; %make multiple of 4
bits = floor(2*rand(N,1));

%number of bits / symbol
%1 or 2 or 4
k = 1;

Ns = N/k;

% en = 0;
% for i=0:15
%     sym_c = desym2complexsym(i,k);
%     en = en + abs(sym_c)^2;
% end
% en/16
snrvec = linspace(10,0,20);
ber = zeros(size(snrvec));
parfor idsnr = 1:length(snrvec)
    snr = snrvec(idsnr);
    rx_data = zeros(N,1);
    for idx = 1:k:N-k+1
        %pick k bits
        tx_bits = bits(idx:idx+k-1,1);
        %interpret as gray code and convert to a symbol number
        sym_de = gray2de(tx_bits,k);
        %convert symbol number to complex sample using constellation
        tx_sym = desym2complexsym(sym_de,k);
        %add noise to sample based on energy in 1 bit and snr 
        rx_sym = add_noise(tx_sym,snr,k);
        %     plot(real(rx_sym),imag(rx_sym),'*')
        %decode received symbol to closest symbol number 
        decoded_sym_de = decode_complex2de(rx_sym,k);
        %convert decoded symbol number to data bits using gray code
        rx_bits = desym2data(decoded_sym_de,k)';
        rx_data(idx:idx+k-1,1) = rx_bits;
        %     hold on
    end
    % bits(:)
    %compare decoded bits to transmitted bits to calculate BER
    ber(idsnr) = 1-sum(rx_data == bits)/N;
end
semilogy(snrvec,ber,'-o','lineWidth',2)
grid on
hold on
semilogy(snrvec,qfunc(sqrt(2*10.^(snrvec/10))),'-x','lineWidth',2)
xlabel("E_b/N_0 dB")
ylabel("P_b")



%%


bits = floor(2*rand(N,1));

%number of bits / symbol
%1 or 2 or 4
k = 2;

Ns = N/k;

% en = 0;
% for i=0:15
%     sym_c = desym2complexsym(i,k);
%     en = en + abs(sym_c)^2;
% end
% en/16
snrvec = linspace(10,0,20);
ber = zeros(size(snrvec));
parfor idsnr = 1:length(snrvec)
    snr = snrvec(idsnr);
    rx_data = zeros(N,1);
    for idx = 1:k:N-k+1
        tx_bits = bits(idx:idx+k-1,1);
        sym_de = gray2de(tx_bits,k);
        tx_sym = desym2complexsym(sym_de,k);
        rx_sym = add_noise(tx_sym,snr,k);
        %     plot(real(rx_sym),imag(rx_sym),'*')
        decoded_sym_de = decode_complex2de(rx_sym,k);
        rx_bits = desym2data(decoded_sym_de,k)';
        rx_data(idx:idx+k-1,1) = rx_bits;
        %     hold on
    end
    % bits(:)
    ber(idsnr) = 1-sum(rx_data == bits)/N;
end
semilogy(snrvec,ber,'-o','lineWidth',2)
grid on
hold on
semilogy(snrvec,qfunc(sqrt(2*10.^(snrvec/10))),'-x','lineWidth',2)
xlabel("E_b/N_0 dB")
ylabel("P_b")



%%


%number of bits in a packet

bits = floor(2*rand(N,1));

%number of bits / symbol
%1 or 2 or 4
k = 4;

Ns = N/k;

% en = 0;
% for i=0:15
%     sym_c = desym2complexsym(i,k);
%     en = en + abs(sym_c)^2;
% end
% en/16
snrvec = linspace(10,0,20);
ber = zeros(size(snrvec));
parfor idsnr = 1:length(snrvec)
    snr = snrvec(idsnr);
    rx_data = zeros(N,1);
    for idx = 1:k:N-k+1
        tx_bits = bits(idx:idx+k-1,1);
        sym_de = gray2de(tx_bits,k);
        tx_sym = desym2complexsym(sym_de,k);
        rx_sym = add_noise(tx_sym,snr,k);
        %     plot(real(rx_sym),imag(rx_sym),'*')
        decoded_sym_de = decode_complex2de(rx_sym,k);
        rx_bits = desym2data(decoded_sym_de,k)';
        rx_data(idx:idx+k-1,1) = rx_bits;
        %     hold on
    end
    % bits(:)
    ber(idsnr) = 1-sum(rx_data == bits)/N;
end
semilogy(snrvec,ber,'-o','lineWidth',2)
grid on
hold on

semilogy(snrvec,1.5/k*erfc(sqrt(k*10.^(snrvec/10)/10)),'-x','lineWidth',2)
legend("BPSK-sim","BPSK-theor","QPSK-sim","QPSK-theor","16-QAM-sim","16-QAM-theor");
xlabel("E_b/N_0 dB")
ylabel("P_b")