cov = zeros(2,2)
Nr = 2;                

for i=1:1000
H = randn(Nr,1)+1j*randn(Nr,1);
H = H./sqrt(2);
%add spatial correlation of 0.8 
H = chol([1,0.8;0.8,1])'*H;
cov = cov + 1/1000*(H*H');
end
cov
