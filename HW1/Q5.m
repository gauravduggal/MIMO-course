close all
clear all
clc

%Receive antennas
Nr = 4;
%Transmit antennas
Nt = 1;
%frequency
f = 5e9;
%wavelength
l = 3e8/f;
%Antenna spacing at rx
d1 = l/2;
%Antenna spacing at tx
d2 = l/2;
%number of multipath
N = 5000;
%sampling frequency
fs = 1e3;
%Time duration
T = 1;
%time axis
t = 0:1/fs:T-1/fs;
%max doppler
fd = 100;
% RIcian K factor
K = 0;
u = 0;
sigma = 90;


T = 0.1;
t = 0:1/fs:T-1/fs;
Nt = 1;
Nr = 2;

vec_d = 0:l/50:5*l-l/50;
c1 = zeros(1,length(vec_d));
for idx = 1:length(vec_d)
    d2 = vec_d(idx);
    H = channel_gaussian_angle_K(Nt, Nr, t, fd, d1,d2,N,l,K,u,sigma);
%     H = channel_uniform_angle_K(Nt, Nr, t, fd, d1,d2,N,l,K);
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
    c1(1,idx) = abs(temp(1,2));%real(cmat(1,101));

end
figure
plot((1:length(vec_d))/50,c1)
ylim([-0.1,1])
xlabel('spacing (wavelengths)')
ylabel('correlation')
title("Spatial correlation vs antenna spacing")
hold on
x= [0.32,1.88,2];
y = [0.8,0.1,0.02];
dy = 0.05;
y_str = string(y);
for id = 1:length(x)
    plot(x(id),y(id),"x",Color="red")
    str = string("corr:"+string(y(id))+", ("+string(x(id))+", "+string(y(id))+")");
    text(x(id),y(id),'\leftarrow '+str,"HorizontalAlignment","left","Color","red","VerticalAlignment","middle",FontSize=14)
end



%%

%Receive antennas
Nr = 2;
%Transmit antennas
Nt = 2;
d2 = 1.88*l;
d1 = 2*l;
%Time duration
T = 0.055;
%time axis
t = 0:1/fs:T-1/fs;
Niter = 10000;
N = 2000;
C1 = zeros(1,Niter);
e1 = zeros(1,Niter);

parfor idt = 1:Niter
    H = channel_gaussian_angle_K(Nt, Nr, t, fd, d1,d2,N,l,K,u,sigma);
    temp = squeeze(H(50,:,:));
    eigvalues = eig(temp'*temp);
%     C1(ctr) = log2(1+10*eigvalues(1)^2)+log2(1+10*eigvalues(2)^2);
    C1(idt) = capacity(temp,10,Nt,Nr);
    e1(1,idt) = abs(max(eigvalues));
end



%tx side 0.8 rx side ~ 0
d2 = 2*l;
d1 = 0.32*l;
C2 = zeros(1,Niter);
e2 = zeros(1,Niter);
parfor idt = 1:Niter
    H = channel_gaussian_angle_K(Nt, Nr, t, fd, d1,d2,N,l,K,u,sigma);
    temp = squeeze(H(50,:,:));
    eigvalues = eig(temp'*temp);

%     C1(ctr) = log2(1+10*eigvalues(1)^2)+log2(1+10*eigvalues(2)^2);
    C2(idt) = capacity(temp,10,Nt,Nr);
    e2(1,idt) = abs(max(eigvalues));
end


d2 = 0.32*l;
d1 = 0.32*l;
%Time duration





C3 = zeros(1,Niter);
e3 = zeros(1,Niter);

parfor idt = 1:Niter
    H = channel_gaussian_angle_K(Nt, Nr, t, fd, d1,d2,N,l,K,u,sigma);
    temp = squeeze(H(50,:,:));
    eigvalues = eig(temp'*temp);
%     C1(ctr) = log2(1+10*eigvalues(1)^2)+log2(1+10*eigvalues(2)^2);
    C3(idt) = capacity(temp,10,Nt,Nr);
    e3(1,idt) = abs(max(eigvalues));
end

d2 = 0.5*l;
d1 = 0.5*l;
Nt = 1;
Nr = 2;
%Time duration





C4 = zeros(1,Niter);
e4 = zeros(1,Niter);

parfor idt = 1:Niter
    H = channel_gaussian_angle_K(Nt, Nr, t, fd, d1,d2,N,l,K,u,sigma);
    temp = squeeze(H(50,:,:));
    eigvalues = eig(temp'*temp);
%     C1(ctr) = log2(1+10*eigvalues(1)^2)+log2(1+10*eigvalues(2)^2);
    C4(idt) = capacity(temp,10,Nt,Nr);
    e4(1,idt) = abs(max(eigvalues));
end



figure;
histogram(e1,100,EdgeColor="none",Normalization='pdf',FaceAlpha=0.9)
hold on
histogram(e2,100,EdgeColor="none",Normalization='pdf',FaceAlpha=0.9)
histogram(e3,100,EdgeColor="none",Normalization='pdf',FaceAlpha=0.9)
xlim([0,10])
legend("tx-0.1, rx-~0","tx-0.8, rx-~0","tx-0.8, rx-0.8")
xlabel("Magnitude of Eigenvalues")
ylabel("probability density")


figure;
hold on
histogram(C3,100,EdgeColor="none",Normalization='cdf',FaceAlpha=0.7)
histogram(C2,100,EdgeColor="none",Normalization='cdf',FaceAlpha=0.7)
histogram(C1,100,EdgeColor="none",Normalization='cdf',FaceAlpha=0.7)
histogram(C4,100,EdgeColor="none",Normalization='cdf',FaceAlpha=0.7)
legend("txcorr-0.1, rxcorr~0, MIMO 2X2","txcorr-0.8, rxcorr~0 ,MIMO 2X2","txcorr-0.8, rxcorr-0.8,MIMO 2X2","txcorr~0, rxcorr~0,SIMO 2X1")
title("Capacity Plots")
xlabel("Capacity (bps/s/Hz")
ylabel("Probability")
ylim([0,1])
mean(C1)
mean(C2)
mean(C3)
mean(C4)


