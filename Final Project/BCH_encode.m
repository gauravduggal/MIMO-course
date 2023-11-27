function coded_bits = BCH_encode(inp,M,k)

    n = 2^M-1;   % Codeword length
%     t = bchnumerr(n,k);
    input_reshaped = gf(reshape(inp,[],k));
    enc = bchenc(input_reshaped,n,k);
    coded_bits_gf = reshape(enc,1,[]);
    coded_bits = double(coded_bits_gf==1);
end