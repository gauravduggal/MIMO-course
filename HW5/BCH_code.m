close all
clear all
clc
M = 3;
n = 2^M-1;   % Codeword length
k = 4;       % Message length
nwords = 50; % Number of words to encode

msgTx = gf(randi([0 1],nwords,k));

%Find the error-correction capability.
t = bchnumerr(n,k);

%Encode the message.
enc = bchenc(msgTx,n,k);

%Corrupt up to t bits in each codeword.
noisycode = enc + randerr(nwords,n,1:t);

%Decode the noisy code.
msgRx = bchdec(noisycode,n,k);

%Validate that the message was properly decoded.
isequal(msgTx,msgRx)

% Increase the number of possible errors, and generate another noisy codeword.
t2 = t + 1;
noisycode2 = enc + randerr(nwords,n,1:t2);

% Decode the new received codeword.
[msgRx2,numerr] = bchdec(noisycode2,n,k);

% Determine if the message was properly decoded by examining the number of corrected errors, numerr. Entries of -1 correspond to decoding failures, which occur when the codeword has more errors than can be corrected for the specified [n,k] pair.
% numerr
% Two of the ten transmitted codewords were not correctly received.
