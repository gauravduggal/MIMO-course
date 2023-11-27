clear all
clc
close all


Ntrials = 100000;

snrvec = 20:-0.5:0;
antvec = [2,4];
kvec = [2,4];
mod_order = ["QPSK","16-QAM"];
for antidx = 1:length(antvec)
    for kidx=1:length(kvec)
        k = kvec(kidx);
        ber =  zeros(size(snrvec));
        parfor snridx=1:length(snrvec)
            tic
            snr = snrvec(snridx);
            ber_snr = 0;
            Nt = antvec(antidx);
            Nr = antvec(antidx);
            MLE_LUT = create_LUT(k,Nt);
                
            for trials=1:Ntrials
                
                H = randn(Nr,Nt)+1j*randn(Nr,Nt);
                H = H./sqrt(2);


                %generate N data bits
                Nbits = Nt*k;
                data_bits = floor(2*rand(Nbits,1));
                bits = data_bits;
                temp = reshape(bits,[k,Nt]);
                sym_dec = bit2int(temp,k)';
                sym = qammod(sym_dec,2^k,'UnitAveragePower',true);
                rx = add_noise(H*sym,snr,k);
%                 rx = H*sym;

                %magic sauce section
                s_est = decode_MLE(rx,H,MLE_LUT);
                s_est_dec = qamdemod(s_est,2^k,'UnitAveragePower',true);
                
                ber_snr = ber_snr +  1/Ntrials*(1 - 1/(k*Nt)*sum(sum(int2bit(s_est_dec,k)==int2bit(sym_dec',k))));

            end
            ber(snridx) = ber_snr;
            toc
        end
        legend_str = string(mod_order(kidx))+" "+string(antvec(antidx))+"x"+string(antvec(antidx))+"MLE-sim";
        semilogy(snrvec,ber,"DisplayName",legend_str);
        hold on
    end
end
ylim([1e-5,1])
grid on
legend
xlabel("E_b/N_0 (dB)")
ylabel("BER")