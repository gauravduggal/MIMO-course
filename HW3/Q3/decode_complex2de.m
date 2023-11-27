function [decoded_sym_de] = decode_complex2de(rx_sym,k)
LUT_16QAM = [-3+3j;
    -1+3j;
    1+3j;
    3+3j;

    3+1j;
    1+1j;
    -1+1j;
    -3+1j;

    -3-1j;
    -1-1j;
    1-1j;
    3-1j;

    3-3j;
    1-3j;
    -1-3j;
    -3-3j]/sqrt(10);

LUT_QPSK = [-1+1j;
    1+1j;
    1-1j;
    -1-1j]/sqrt(2);

LUT_BPSK = [1;
    -1]/sqrt(2);

min_d = inf;
decoded_sym_de = 0;
switch(k)
    case 4
        for i=1:2^k
            d = abs(rx_sym - LUT_16QAM(i,1))^2;
            if min_d>d
                min_d = d;
                decoded_sym_de = i-1;
            end
        end
    case 2
        for i=1:2^k
            d = abs(rx_sym - LUT_QPSK(i,1))^2;
            if min_d>d
                min_d = d;
                decoded_sym_de = i-1;
            end
        end
    case 1
        for i=1:2^k
            d = abs(rx_sym - LUT_BPSK(i,1))^2;
            if min_d>d
                min_d = d;
                decoded_sym_de = i-1;
            end
        end
    otherwise
        disp("k = 1,2 or 4")
end
end