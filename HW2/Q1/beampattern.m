function [bp] = beampattern(w,N,lambda,d)
thetavec = 0:360;
bp = zeros(size(thetavec));
for idx = 1:length(thetavec)
    bp(idx) = abs(w'/norm(w)*array_response(thetavec(idx),N,lambda,d));
end
figure;
polarpattern(thetavec,bp,'AngleAtTop',0,'AngleLim', [-180 180]);

end

