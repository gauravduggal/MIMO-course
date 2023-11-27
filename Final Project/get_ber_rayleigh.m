function [pb] = get_ber_rayleigh(snr_linear,klist)
pb = zeros(length(snr_linear),1);

for snridx = 1:length(snr_linear)
    k = klist(snridx);
    eb_no = 10*log10(1/k*snr_linear(snridx));
    switch k
        case 0
            pb(snridx,1) = 0;
        case 1
            pb(snridx,1) = berfading(eb_no,'psk',2,1);
        case 2
            pb(snridx,1) = berfading(eb_no,'qam',4,1);
        case 3
            pb(snridx,1) = berfading(eb_no,'qam',8,1);            
        case 4
            pb(snridx,1) = berfading(eb_no,'qam',16,1);
        case 5
            pb(snridx,1) = berfading(eb_no,'qam',32,1);            
        case 6
            pb(snridx,1) = berfading(eb_no,'qam',64,1);
        case 8
            pb(snridx,1) = berfading(eb_no,'qam',256,1);
        
        otherwise
            pb = 0;
    end
end
end