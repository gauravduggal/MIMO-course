function [de] = gray2de(bits,k)
% k = length(bits);
de = -1;
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
        for i=1:2^k
            gray = LUT_16QAM(i,:);
            if (sum(bits==gray')==k)
                de = i-1;
                return
            end
        end
    case 2
        for i=1:2^k
            gray = LUT_QPSK(i,:);
            if (sum(bits==gray')==k)
                de = i-1;
                return
            end
        end
    case 1
        for i=1:2^k
            gray = LUT_BPSK(i,:);
            if (sum(bits==gray')==k)
                de = i-1;
                return
            end
        end
    otherwise
        disp("k = 1,2 or 4")
end

end
