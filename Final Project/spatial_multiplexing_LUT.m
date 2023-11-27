close all
clear all
clc


Nt = 2;
Nr = 2;
Ntrials = 1000;
% kvec = [1,2,3,4,5,6];
kvec = 2*ones(1,6);
snrvec = 30:-2:0;
ber = zeros(size(snrvec));
for snridx = 1:length(snrvec)
    snr = 10^(snrvec(snridx)/10);
    for trial = 1:Ntrials

        k_tx1 = kvec(floor(6*rand(1)+1));
        bits1 = floor(2*rand(k_tx1,1));
        s1_dec = bit2int(bits1,k_tx1);
        if k_tx1==1
            s1 = bits1;
        else
            s1 = qammod(s1_dec,2^k_tx1,'UnitAveragePower',true);
        end



        k_tx2 = kvec(floor(6*rand(1)+1));
        bits2 = floor(2*rand(k_tx2,1));
        s2_dec = bit2int(bits2,k_tx2);
        if k_tx2==1
            s2 = bits2;
        else
            s2 = qammod(s2_dec,2^k_tx2,'UnitAveragePower',true);
        end
        bits_in = [bits1;bits2];


        h11 = randn(1) + 1j*randn(1);
        h11 = h11/sqrt(2);
        h12 = randn(1) + 1j*randn(1);
        h12 = h12/sqrt(2);

        h21 = randn(1) + 1j*randn(1);
        h21 = h21/sqrt(2);
        h22 = randn(1) + 1j*randn(1);
        h22 = h22/sqrt(2);
        H = [h11,h12;h21, h22];

        s = [s1;s2];
        n = (randn(2,1)+1j*randn(2,1))/sqrt(2*snr);
        r = H*s+n;
        z = (H'*H+1/snr*eye(2))^-1*H'*r;

        rxbits1 = int2bit(qamdemod(z(1),2^k_tx1,'UnitAveragePower',true),k_tx1);

        rxbits2 = int2bit(qamdemod(z(2),2^k_tx2,'UnitAveragePower',true),k_tx2);

        bits_out = [rxbits1;rxbits2];

        ber(snridx) = ber(snridx) + 1/Ntrials*mean(bits_in~=bits_out);


    end
end
 semilogy(snrvec,ber)
 hold on
 semilogy(snrvec,get_ber_rayleigh(10.^(snrvec/10),6*ones(length(snrvec),1)))