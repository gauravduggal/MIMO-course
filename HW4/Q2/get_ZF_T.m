function [T] = get_ZF_T(H)
T = (H'*H)^(-1)*H';
end