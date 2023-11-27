clear all
clc
close all


Ntrials = 10000;
snrvec = 20:-0.5:0;
antvec = [2,4];
kvec = [1,2,4];
mod_order = ["BPSK","QPSK","16-QAM"];
for antidx = 1:length(antvec)
    for kidx=1:length(kvec)
        k = kvec(kidx);
        ber =  zeros(size(snrvec));
        parfor snridx=1:length(snrvec)
            tic
            snr = snrvec(snridx);
            ber_snr = 0;
            for trials=1:Ntrials
                Nt = antvec(antidx);
                Nr = antvec(antidx);

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
%                             rx = H*sym;

                %Plain Zero Forcing
                T = get_ZF_T(H);
                z = T*rx;
                s_est = qamdemod(z,2^k,'UnitAveragePower',true);
                ber_snr = ber_snr +  1/Ntrials*(1 - 1/(k*Nt)*sum(sum(int2bit(s_est,k)==int2bit(sym_dec,k))));
%                 s_est
%                 sym_dec'
%                 ber_snr = ber_snr +  1/Ntrials*(1 - 1/(k*Nt)*sum(sum(int2bit(sort(s_est),k)==int2bit(sort(sym_dec)',k))));
%                             s_est
%                             sym_dec(decoding_order)'
                %             decoding_order


            end
            ber(snridx) = ber_snr;
            toc
        end
        legend_str = string(mod_order(kidx))+" "+string(antvec(antidx))+"x"+string(antvec(antidx))+"ZF";
        semilogy(snrvec,ber,'-x',"DisplayName",legend_str);
        hold on
    end
end
ylim([1e-5,1])
grid on
xlabel("E_b/N_0")
ylabel("BER")
semilogy(snrvec,BER_ZF_BPSK(2,snrvec),"-o","DisplayName","BPSK-2x2,ZF-theoretical");
hold on
semilogy(snrvec,BER_ZF_BPSK(4,snrvec),"-^","DisplayName","BPSK-4x4-ZF-theoretical");
legend
