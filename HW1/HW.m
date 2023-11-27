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
d1 = l/2;
%Antenna spacing at tx
d2 = l/2;
%number of multipath
N = 200;
%sampling frequency
fs = 1e3;
%Time duration
T = 0.1;
%time axis
t = 0:1/fs:T-1/fs;
%max doppler
fd = 100;

%% Q1 a).
Nt = 1;
Nr = 2;
d1 = 0;
d2 = l/2;
fd = 100;


vec_d = 0:l/50:3*l-l/50;
c = zeros(1,length(vec_d));
for idx = 1:length(vec_d)
    d2 = vec_d(idx);
    H = channel_uniform_angle(Nt, Nr, t, fd,d1,d2,N,l);
%     H = channel_gaussian_angle(Nt, Nr, t, fd, u, sigma,d1,d2,N,l);
    H = reshape(H(:,1,:),[Nr,length(t)]);

    c(1,idx) = corr(abs(H(1,:))',abs(H(2,:))');
end
plot((1:length(vec_d))/50,c)
xlabel('spacing (wavelengths)')
ylabel('correlation')

H = channel_uniform_angle(Nt, Nr, t, fd,d1,d2,N,l);
H = reshape(H(:,1,:),[Nr,length(t)]);
corr(abs(H(1,:))',abs(H(2,:))')
semilogy(t,abs(H(1,:)))
hold on
semilogy(t,abs(H(2,:)))
xlabel('time (s)')
ylabel('signal amplitude (db)')


%% Q1 b).
Nt = 1;
Nr = 2;
d1 = l/2;
d2 = l/2;
fd = 100;
T = 10;
t = 0:1/fs:T-1/fs;
H = channel_uniform_angle(Nt, Nr, t, fd,d1,d2,N,l);
H = reshape(H(:,1,:),[Nr,length(t)]);
corr(abs(H(1,:))',abs(H(2,:))')
sum(abs(H(1,:)).^2)/T
sum(abs(H(2,:)).^2)/T

%% Q1 c). Temporal correlation
Nt = 1;
Nr = 1;
d1 = l/2;
d2 = l/2;
fd = 100;
T = 30;
t = 0:1/fs:T-1/fs;
H = channel_uniform_angle(Nt, Nr, t, fd,d1,d2,N,l);
H = reshape(H(:,1,:),[Nr,length(t)]);

hrx1 = abs(H);
samples_axis = 1:T*fs;
signal_length_samples = 0.1*fs;

for tau = 1:signal_length_samples
%     tau
    ctr = 0;
    c = 0;
    for sample=1:signal_length_samples:length(samples_axis)-2*signal_length_samples
        rt1 = hrx1(sample:sample+signal_length_samples-1);
%         sample+tau+signal_length_samples-1
        rt2 = hrx1(sample+tau:sample+tau+signal_length_samples-1);
        temp = rt1*rt2'/(sqrt(sum(abs(rt1).^2))*sqrt(sum(abs(rt2).^2)));
%         temp = corr(rt1',rt2');

        ctr = ctr + 1;
        c = c + temp;
    end
    c = c / ctr;
    plot(tau,c,"*")
    hold on
end
ylim([-1,1])

%% Q 1d). fd is changed
Nt = 1;
Nr = 1;
d1 = l/2;
d2 = l/2;
fd = 15;
T = 30;
t = 0:1/fs:T-1/fs;
H = channel_uniform_angle(Nt, Nr, t, fd,d1,d2,N,l);
H = reshape(H(:,1,:),[Nr,length(t)]);

hrx1 = abs(H);
samples_axis = 1:T*fs;
signal_length_samples = 0.1*fs;

for tau = 1:signal_length_samples
%     tau
    ctr = 0;
    c = 0;
    for sample=1:signal_length_samples:length(samples_axis)-2*signal_length_samples
        rt1 = hrx1(sample:sample+signal_length_samples-1);
%         sample+tau+signal_length_samples-1
        rt2 = hrx1(sample+tau:sample+tau+signal_length_samples-1);
        temp = rt1*rt2'/(sqrt(sum(abs(rt1).^2))*sqrt(sum(abs(rt2).^2)));
%         temp = corr(rt1',rt2');

        ctr = ctr + 1;
        c = c + temp;
    end
    c = c / ctr;
    plot(tau,c,"*")
    hold on
end
ylim([-1,1])



%% Q 2 Spatial correlation

Nt = 1;
Nr = 2;
d1 = 0;
d2 = l/2;

%Time duration
T = 100;
%time axis
t = 0:1/fs:T-1/fs;
%max doppler
fd = 100;



H = channel_uniform_angle(Nt, Nr, t, fd,d1,d2,N,l);
H = reshape(H(:,1,:),[Nr,length(t)]);
c = [];
for sample = 1:length(t)-signal_length_samples
c = [c, corr(abs(H(1,sample:sample+signal_length_samples))',abs(H(2,sample:sample+signal_length_samples))')];
end
plot(t(1:end-signal_length_samples),c)

%% Q3 

%% Q4 a).
Nt = 1;
Nr = 2;
N = 100;
d1 = 0;
d2 = l/2;
fd = 100;
H = channel_uniform_angle(Nt, Nr, t, fd,d1,d2,N,l);
H = reshape(H(:,1,:),[Nr,length(t)]);

rt = abs(H(1,:));

histogram(rt,100)

%% Q4 b). 
Nt = 1;
Nr = 2;
N = 100;
d1 = 0;
d2 = l/2;
fd = 100;
K = 0;
H = channel_uniform_angle_K(Nt, Nr, t, fd,d1,d2,N,l,K);
H = reshape(H(:,1,:),[Nr,length(t)]);

rt = abs(H(1,:));

histogram(rt,100)

