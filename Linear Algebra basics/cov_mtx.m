function [C]  = cov_mtx(s)
%nomalises covariance matrix so all variances are 1
    %covariance matx
    temp = s'*s;
    %making variance = 1
    std_dev = (diag(temp)).^-0.5;
    C = diag(std_dev)*temp*diag(std_dev);
end