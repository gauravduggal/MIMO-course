function [k_list, avg_ber,Ne,enabled_scs,snr_db_w_null] =  get_bitloading(snr_linear, target_ber)
N = length(snr_linear);
pow_list= ones(N,1);
%assuming tx power on each subcarrier is 1
noise_pow_list = 1./snr_linear;

ber_per_sc = get_ber_awgn2(snr_linear,6*ones(N,1));
k_list = 6*ones(N,1);
enabled_scs = ones(N,1);
Ne = sum(enabled_scs);
avg_ber = sum(ber_per_sc.*k_list)/(sum(k_list));
while(avg_ber> target_ber)
    
    [~,idx] = max(ber_per_sc);
    k_list(idx) = get_lower_modulation_scheme(k_list(idx));
%     if k_list(idx)==0
%         enabled_scs = k_list>0;
%         Ne = sum(enabled_scs);
%         pow_list(enabled_scs) = ones(Ne,1)*N/Ne;
%         pow_list(~enabled_scs) = 0;
%         snr_linear = pow_list./noise_pow_list;
%     end
    enabled_scs = k_list>0;
    Ne = sum(enabled_scs);
    ber_per_sc = get_ber_awgn2(snr_linear,k_list);
    avg_ber = sum(ber_per_sc.*k_list)/(sum(k_list)); 
    snr_db_w_null = 10*log10(snr_linear);
end