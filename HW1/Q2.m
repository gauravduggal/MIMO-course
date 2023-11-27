close all
clear all
clc


%Receive antennas
Nr = 2;
%Transmit antennas
Nt = 1;
%frequency
f = 5e9;
%wavelength
l = 3e8/f;
%Antenna spacing at rx
d1 = l/100;
%Antenna spacing at tx
d2 = l/4;
%number of multipath
N = 1000;
%sampling frequency
fs = 1e3;
%Time duration
T = 0.1;
%time axis
t = 0:1/fs:T-1/fs;
%max doppler
fd = 100;
%Rician K factor
K = 0;


% Q2
% c = zeros(1,length(t));
% H = channel_uniform_angle_K(Nt, Nr, t, fd,d1,d2,N,l,K);
% for sample=1:1:length(t)-1000
%     temp=corrcoef(H(sample:sample+100,:));
%     c(1,sample) = temp(1,2);
% end
% plot(abs(c(1:length(t)-1000)))
% mean(abs(c))

c = zeros(1,length(t));
Niter = 1000;

%Assuming ergodicity we create a channel 100 times and calculation spatial
%correlation between antennas
parfor idx = 1:Niter
    H = channel_uniform_angle_K(Nt, Nr, t, fd,d1,d2,N,l,K);
    temp=corrcoef(H);
    c(1,idx) = temp(1,2);% H(1,1)*H(1,2)'/2;
end

plot((1:Niter)/10,abs(c))
hold on
d2 = l/2;
c = zeros(1,length(t));

%Assuming ergodicity we create a channel 100 times and calculation spatial
%correlation between antennas
parfor idx = 1:Niter
    H = channel_uniform_angle_K(Nt, Nr, t, fd,d1,d2,N,l,K);
    temp=corrcoef(H);
    c(1,idx) = temp(1,2);% H(1,1)*H(1,2)'/2;
end

plot((1:Niter)/10,abs(c))
ylabel("Spatial correlation")
xlabel("time / expt number")
title("Spatial correlation over 100 seconds")
hold on
legend("d=l/4","d=l/2")

%%
T = 10;
N = 1000;
%time axis
t = 0:1/fs:T-1/fs;

vec_d = l/50:l/50:5*l-l/50;
c1 = zeros(1,length(vec_d));

parfor idx = 1:length(vec_d)
    d2 = vec_d(idx);
%     H = channel_gaussian_angle_K(Nt, Nr, t, fd, d1,d2,N,l,K,u,sigma);
    H = channel_uniform_angle_K(Nt, Nr, t, fd, d1,d2,N,l,K);
    temp = corrcoef(H);
    
% %     H = transpose(H);
%     Hvec = H(:);
%     cmat = Hvec*Hvec';
%     c = 0;
%     for sample=1:length(t)
%         Hinst = H(sample,:,:);
%         temp = Hinst'*Hinst;
%         c = c + temp(1,2);
%     end
%     c = c/length(t);
%     c1(1,idx) = abs(c);
    c1(1,idx) = real(temp(1,2));%real(cmat(1,101));


end
figure
plot((1:length(vec_d))/50,c1)
ylim([-0.5,1])
xlabel('spacing (wavelengths)')
ylabel('correlation')
title("Spatial correlation vs antenna spacing")



hold on
sp = besselj(0,vec_d*2*pi/l);
plot(vec_d/l,sp,'-x');
legend("simulated","theoretical")

