function [decoded_sym] = receive(rx,w, nbits_sym, AoA,lambda,d)
N = length(w);
ar = array_response(AoA,N,lambda,d);
z = w'*ar*rx;
decoded_sym = qamdemod(z,nbits_sym);
end

