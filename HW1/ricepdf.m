function y = ricepdf(x, K,sigma)
%%https://ieeexplore-ieee-org.ezproxy.lib.vt.edu/stamp/stamp.jsp?tp=&arnumber=1210745
%%The Ricean K Factor: Estimation and Performance Analysis

try
    y = (2*(K+1)*x/sigma) .*...
        exp(-K*ones(size(x))-(K+1)/sigma*x.^2) .*...
        besseli(0, 2*x .* sqrt(K*(K+1)/sigma));
        % besseli(0, ...) is the zeroth order modified Bessel function of
        % the first kind. (see help bessel)
    y(x <= 0) = 0;
catch
    error('ricepdf:InputSizeMismatch',...
        'Non-scalar arguments must match in size.');
end
