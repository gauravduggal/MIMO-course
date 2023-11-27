clear all
clc
close all


Ntrials = 100000;
snrvec = 20:-0.5:0;

corrvec = [0,0.7];
mod_order = ["BPSK"];
for corridx=1:length(corrvec)
    corr = corrvec(corridx);
    k = 1;
    Nt = 2;
    Nr = 2;

    ber =  zeros(size(snrvec));
    parfor snridx=1:length(snrvec)
        snr = snrvec(snridx);
        ber_snr = 0;
        for trials=1:Ntrials
            
            H = randn(Nr,Nt)+1j*randn(Nr,Nt);
            H = H./sqrt(2);
            H = chol([1,corr;corr,1])'*H;

            %generate N data bits
            Nbits = Nt*k;
            data_bits = floor(2*rand(Nbits,1));
            bits = data_bits;
            temp = reshape(bits,[k,Nt]);
            sym_dec = bit2int(temp,k)';
            sym = qammod(sym_dec,2^k,'UnitAveragePower',true);
            rx = add_noise(H*sym,snr,k);
            %                             rx = H*sym;

            %ZF
            s_est = zeros(1,Nt);
            sym_to_tx_ant_LUT = 1:Nt;
            decoding_order = zeros(1,Nt);
            ctr = 1;
            %magic sauce section
            while (ctr<=Nt)
                %find row with min noise amplification
                %                     T = get_ZF_T(H);
                T = get_MMSE_T(H,Nt,snr,k);
                norm2 = abs(T).^2;
                norm2 = sum(norm2,2);
                [minv,minidx] = min(norm2);
                z = T(minidx,:)*rx;
                s_hat_dec = qamdemod(z,2^k,'UnitAveragePower',true);
                s_hat = qammod(s_hat_dec,2^k,'UnitAveragePower',true);
                decoding_order(ctr) = sym_to_tx_ant_LUT(minidx);
                %                     sym_to_tx_ant_LUT
                sym_to_tx_ant_LUT(minidx) = [];
                rx = rx - s_hat*H(:,minidx);
                H(:,minidx) = [];%zeros(Nr,1);
                s_est(ctr) = s_hat_dec;
                ctr = ctr + 1;
            end
            %             sort(s_est)
            %             sort(sym_dec')
            ber_snr = ber_snr +  1/Ntrials*(1 - 1/(k*Nt)*sum(sum(int2bit(s_est,k)==int2bit(sym_dec(decoding_order)',k))));
            %                 s_est
            %                 sym_dec'
            %                 ber_snr = ber_snr +  1/Ntrials*(1 - 1/(k*Nt)*sum(sum(int2bit(sort(s_est),k)==int2bit(sort(sym_dec)',k))));
            %                             s_est
            %                             sym_dec(decoding_order)'
            %             decoding_order


        end
        ber(snridx) = ber_snr;
    end
    legend_str = "corr-"+string(corr)+" "+string(mod_order)+" "+string(Nt)+"x"+string(Nr)+"sim";
    semilogy(snrvec,ber,"DisplayName",legend_str);
    hold on
end
ylim([1e-5,1])
grid on
legend