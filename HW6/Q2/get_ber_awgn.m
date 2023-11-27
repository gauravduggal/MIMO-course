function [pb] = get_ber_awgn(snr_db,k)
pb = zeros(1,length(snr_db));
eb_no = 10*log10(1/k*10.^(snr_db/10));
for snridx = 1:length(snr_db)
    switch k
        case 0
            pb(1,snridx) = 0;
        case 1
            pb(1,snridx) = berawgn(eb_no(1,snridx),'psk',2,'nondiff');
        case 2
            pb(1,snridx) = berawgn(eb_no(1,snridx),'qam',4,'nondiff');
        case 4
            pb(1,snridx) = berawgn(eb_no(1,snridx),'qam',16,'nondiff');
        case 6
            pb(1,snridx) = berawgn(eb_no(1,snridx),'qam',64,'nondiff');
        otherwise
            pb = nan;
    end
end
end