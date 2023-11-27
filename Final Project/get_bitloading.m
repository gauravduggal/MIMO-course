function [k_list_overall, avg_ber,Ne,enabled_scs,snr_db_w_null] =  get_bitloading(snr_linear, target_ber,data_scs)
N = length(snr_linear);
Nd = length(data_scs);
% pow_list= ones(N,1);
% %assuming tx power on each subcarrier is 1
% noise_pow_list = 1./snr_linear;
kmax = 6;
ber_per_sc = get_ber_awgn2(snr_linear(data_scs),kmax*ones(Nd,1));
k_list = kmax*ones(Nd,1);
enabled_scs = ones(Nd,1);
Ne = sum(enabled_scs);
avg_ber = sum(ber_per_sc.*k_list)/(sum(k_list));

% max_ber_sc = max(ber_per_sc);
% while(max_ber_sc>target_ber)
%     [~,idx] = max(ber_per_sc);
%     k_list(idx) = get_lower_modulation_scheme(k_list(idx));
%     ber_per_sc(idx) = get_ber_rayleigh(snr_linear(idx),k_list(idx));
%     max_ber_sc = max(ber_per_sc);
%     enabled_scs = k_list>0;
% end
% avg_ber = sum(ber_per_sc.*k_list)/(sum(k_list)); 


while(avg_ber> target_ber)
    
    [~,idx] = max(ber_per_sc);
    k_list(idx) = get_lower_modulation_scheme(k_list(idx));

    enabled_scs = k_list>0;
    ber_per_sc(idx) = get_ber_awgn2(snr_linear(idx),k_list(idx));
    avg_ber = sum(ber_per_sc.*k_list)/(sum(k_list)); 
  
end
    k_list_overall = zeros(N,1);
    k_list_overall(data_scs) = k_list; 
    Ne = sum(enabled_scs);
    snr_db_w_null = 10*log10(snr_linear);

end