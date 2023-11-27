function [data_bits,data_symbols] = rx_get_data_sc(rx1_equalized_noisy,rx2_equalized_noisy, k_per_sc_tx1,k_per_sc_tx2,data_scs)
%load data subcarriers with the right modulation order
%total scs
N = length(k_per_sc_tx1);
% tx1_complex_sym = zeros(N,1);
% tx2_complex_sym = zeros(N,1);
%data
data_bits = [];
data_symbols = zeros(N,2);
for idx = 1:length(data_scs)
    data_scidx = data_scs(idx);

    k_tx1 = k_per_sc_tx1(data_scidx);
    if k_tx1 ==0
        continue
    end
    sym_dec = qamdemod(rx1_equalized_noisy(data_scidx),2^k_tx1,'UnitAveragePower',true);
    bits = int2bit(sym_dec,k_tx1);
    data_bits = [data_bits; bits];
    data_symbols(data_scidx,1) = sym_dec;

    k_tx2 = k_per_sc_tx2(data_scidx);
    if k_tx2 ==0
        continue
    end
    sym_dec = qamdemod(rx2_equalized_noisy(data_scidx),2^k_tx2,'UnitAveragePower',true);
    bits = int2bit(sym_dec,k_tx2);
    data_bits = [data_bits; bits];
    data_symbols(data_scidx,2) = sym_dec;
end

end