function [arr] = array_response(theta,N,lambda,d)
%TX_ARRAY Summary of this function goes here
%   Detailed explanation goes here
arr = exp(1j*2*pi/lambda*d*sind(theta)*(1:N)');
end

