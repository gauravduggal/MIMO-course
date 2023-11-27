function [gray] = desym2data(sym_de,k)
LUT_16QAM = [0	0	0	0;
    0	0	0	1	;
    0	0	1	1	;
    0	0	1	0	;
    0	1	1	0	;
    0	1	1	1	;
    0	1	0	1	;
    0	1	0	0	;
    1	1	0	0	;
    1	1	0	1;
    1	1	1	1;
    1	1	1	0;
    1	0	1	0;
    1	0	1	1;
    1	0	0	1;
    1	0	0	0;];


LUT_QPSK = [0	0;
            0	1	;
            1	1	;
            1	0	;];

LUT_BPSK = [0;
            1];
    
switch(k)
    case 4 
        gray = LUT_16QAM(sym_de+1,:);
    case 2
        gray = LUT_QPSK(sym_de+1,:);
    case 1
        gray = LUT_BPSK(sym_de+1,:);
    otherwise
        disp("k = 1,2 or 4")
end
end

