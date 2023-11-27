function [complex_sym] = desym2complexsym(sym_de,k)
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
            -1];
            
    
switch(k)
    case 4 
        complex_sym = LUT_16QAM(sym_de+1,1);
    case 2
        complex_sym = LUT_QPSK(sym_de+1,1);
    case 1
        complex_sym = LUT_BPSK(sym_de+1,1);
    otherwise
        disp("k = 1,2 or 4")
end


end

