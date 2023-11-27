function [ x,y ] = GolayPair( N )
%GOLAYPAIR Summary of this function goes here
%   http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.1004.741&rep=rep1&type=pdf
x = zeros(1,N);
y = zeros(1,N);
for r=0:N-1
    x(r+1) = hn(0,r,1,N);
    y(r+1) = hn(N/2,r,1,N);
end

end
%% read about k. Currently set to 1
%% This function generates the i,j element of Golay Matrix.
function h=hn(i, j,k,n)
j=bitand(bitxor(bitshift(j,1),bitxor(i,k)),j);
for i=1:n-1
j=bitxor(bitshift(j,-1),bitand(j,1));
end
h=(-1)^j;
end