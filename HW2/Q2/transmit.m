function [tx] = transmit(x,w,AoD,lambda,d)
N = length(w);
ar = array_response(AoD,N,lambda,d);
tx = ar*x;
end

