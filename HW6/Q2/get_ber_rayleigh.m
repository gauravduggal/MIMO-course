function [pb] = get_ber_rayleigh(snr_db,k)
pb = zeros(1,length(snr_db));
eb_no = 10*log10(1/k*10.^(snr_db/10));
for snridx = 1:length(snr_db)
    switch k
        case 1
            pb(1,snridx) = berfading(eb_no(1,snridx),'psk',2,1);
        case 2
            pb(1,snridx) = berfading(eb_no(1,snridx),'qam',4,1);
        case 4
            pb(1,snridx) = berfading(eb_no(1,snridx),'qam',16,1);
        case 6
            pb(1,snridx) = berfading(eb_no(1,snridx),'qam',64,1);
        case 8
            pb(1,snridx) = berfading(eb_no(1,snridx),'qam',256,1);
        
        otherwise
            pb = 0;
    end
end
end