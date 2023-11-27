function [value] = gcna(c,n,a)

u = sqrt(c/(c+2*a));
temp = 0;
for k=0:n-1
    temp = temp + nchoosek(2*k,k)*((1-u^2)/4)^k;
end
value = 1/(2*a^n)*(1-u*temp);
return