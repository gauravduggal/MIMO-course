function decoded_bits = BCH_decode(noisy_inp,M,k)

    n = 2^M-1;   % Codeword length
    input_reshaped = gf(reshape(noisy_inp,[],2^M -1));
    dec = bchdec(input_reshaped,n,k);
    decoded_bits_gf = reshape(dec,1,[]);
    decoded_bits = double(decoded_bits_gf==1);
end