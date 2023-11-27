clear all
clc
close all


Ntrials = 1000;
snrvec = 20:-0.5:0;

k = 2;
Nt = 4;
Nr = 4;
corr = 0;
ber =  zeros(size(snrvec));
avg_n= zeros(size(snrvec));
avg_C= zeros(size(snrvec));
parfor snridx=1:length(snrvec)
    tic
    snr = snrvec(snridx);
    ber_snr = 0;
    
    for trials=1:Ntrials

        H = randn(Nr,Nt)+1j*randn(Nr,Nt);
        H = H./sqrt(2);
        corr_matrix =  0.99*~eye(Nr);
        corr_matrix = corr_matrix + eye(Nr);
    
        H = chol([corr_matrix])'*H;

        %generate N data bits
        Nbits = Nt*k;
        data_bits = floor(2*rand(Nbits,1));
        bits = data_bits;
        temp = reshape(bits,[k,Nt]);
        sym_dec = bit2int(temp,k)';
        sym = qammod(sym_dec,2^k,'UnitAveragePower',true);


        %waterfilling magic
        r = rank(H);
        flag =true;
        p = 1;
        while(flag==true)
            [U,S,V] = svd(H);
            %get eigenvalues from singular values
            lambdavec = diag(S).^2;
            %waterlevel
            u = Nt/(r-p+1)*(1+1/(k*10^(snr/10))*sum(1./lambdavec(1:(r-p+1))));
            %gain per mode
            n = u-1/(k*10^(snr/10)) * Nt./lambdavec((1:(r-p+1)));
            if (sum(n<0)==0)
                flag=false;
            else
                p = p + 1;
            end
        end
        temp = zeros(1,Nt);
        N_modes = length(n);
        temp(1:length(n)) = sqrt(n);
        T = diag(temp);
        rx = add_noise(H*V*T*sym,snr,k);
%         rx = H*V*T*sym;
%         sum(abs(rx).^2)
%         
%         %magic sauce section
        z = (U'*rx)./(diag(S));
        z = z(1:N_modes);
        s_est = qamdemod(z,2^k,'UnitAveragePower',true)';
        avg_n(snridx) = avg_n(snridx) + N_modes/(Ntrials);
        ber_snr = ber_snr +  1/Ntrials*(1 - 1/(k*N_modes)*sum(sum(int2bit(s_est(1:N_modes),k)==int2bit(sym_dec(1:N_modes)',k))));
        for nidx = 1:N_modes
            avg_C(snridx) = avg_C(snridx) + (1/Ntrials)* log2(1+k*10^(snr/10)*n(nidx)/Nt*S(nidx,nidx)^2)  
        end
    end
    ber(snridx) = ber_snr;
    toc
end
legend_str = string(Nt)+"x"+string(Nr)+"sim";
semilogy(snrvec,ber,"DisplayName",legend_str);
xlabel("E_b/N_0")
ylabel("BER")
ylim([1e-5,1])
grid on
legend
figure;
plot(snrvec,avg_n)
grid on
xlabel("E_b/N_0")
ylabel("Average spatial modes")
ylim([0,Nt+1])
figure;
plot(10*log10(k*10.^(snrvec/10)),avg_C)
grid on
xlabel("E_b/N_0")
ylabel("Average Capacity over trials")