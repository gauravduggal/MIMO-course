function [T] = get_MMSE_T(H,Nt,eb_n0_db,k)
eb_n0 = 10^(eb_n0_db/10);
T = (H'*H+(Nt*eye(size(H,2))/(k*eb_n0)))^(-1)*H';
end