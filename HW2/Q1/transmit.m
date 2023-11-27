function [tx] = transmit(x,w,AoD,lambda,d)
N = length(w);
ar = array_response(AoD,N,lambda,d);
tx = w'*ar*x;
end

