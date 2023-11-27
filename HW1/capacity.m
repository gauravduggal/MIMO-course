function [c] = capacity(H,sigma,Nt,Nr)
%CAPACITY Summary of this function goes here
%   Detailed explanation goes here
% c = log2(det(eye(Nr)+gamma*H*H'/Nt));
c = log2(abs(det(eye(Nr)  + (sigma * ((H *H') / Nt)) )));
end

