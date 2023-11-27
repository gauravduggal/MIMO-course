function [s_hat] = decode_MLE(r,H,LUT)
Nt = size(H,2);
symbol_amplitude = 1;%sqrt(1/Nt);
N_sym = size(LUT,1);
d = zeros(N_sym,1);
for i=1:N_sym
    d(i) = sum(abs(r-symbol_amplitude*H*transpose(LUT(i,:))).^2);
end
[~,minidx] = min(d);
s_hat = LUT(minidx,:);

end