function [Nbits] = get_throughput_rayleigh(snr,N,ofdm_symbols,k,target_ber)
    ber = get_ber_rayleigh(snr,k);
    ber(ber>target_ber) = 1; 
%     for beridx = 1:length(ber)
%         for nsc = N:-1:1
%             if ber(beridx)>target_ber
%                 if (ber(beridx)*nsc/N < target_ber)
%                     ber(beridx) = ber(beridx)*nsc/N;
%                     break;
%                 end
%             end
%         end
%     end
%     ber
    Nbits = (1-ber)*ofdm_symbols*N*k;

end