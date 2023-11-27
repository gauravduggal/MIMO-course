function [c] = corrvec(x,y)
%CORRVEC Summary of this function goes here
%   Detailed explanation goes here
if length(x) ~= length(y)
    disp("x and y need to be same size")
    return
end
lenx = sqrt(sum(x.^2));
leny = sqrt(sum(y.^2));
c = abs(transpose(x)*conj(y)/(lenx*leny));

end

