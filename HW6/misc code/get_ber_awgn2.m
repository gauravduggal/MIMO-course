function [pb] = get_ber_awgn2(snr_linear,klist)
pb = zeros(length(snr_linear),1);

for snridx = 1:length(snr_linear)
    k = klist(snridx);
    eb_no = 10*log10(1/k*snr_linear(snridx));
    switch k
        case 0
            pb(snridx,1) = 0;
        case 1
            pb(snridx,1) = berawgn(eb_no,'psk',2,'nondiff');
        case 2
            pb(snridx,1) = berawgn(eb_no,'qam',4,'nondiff');
        case 4
            pb(snridx,1) = berawgn(eb_no,'qam',16,'nondiff');
        case 6
            pb(snridx,1) = berawgn(eb_no,'qam',64,'nondiff');
        otherwise
            pb = nan;
    end
end
end