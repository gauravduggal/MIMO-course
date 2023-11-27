function [LUT] = create_LUT(k,Nt)
M = 2^k;
LUT = zeros(M^Nt,Nt);
for i=1:M^Nt
   temp = i-1;
   for dim=1:Nt
        LUT(i,dim) = mod(temp,M);
        temp = floor(temp/M);
   end
end

LUT = qammod(LUT,M,'UnitAveragePower',true);
end