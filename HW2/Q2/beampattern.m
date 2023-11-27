function [bp,thetavec] = beampattern(w,N,lambda,d,fig)
thetavec = -180:180;
bp = zeros(size(thetavec));
for idx = 1:length(thetavec)
    bp(idx) = abs(w'/norm(w)*array_response(thetavec(idx),N,lambda,d));
end

% plot(thetavec,bp);
if fig==true
    figure;
    polarpattern(thetavec,bp,'AngleAtTop',0,'AngleLim', [-180 180],'AngleDirection', 'cw');
end
% polarpattern(thetavec,bp);
end

