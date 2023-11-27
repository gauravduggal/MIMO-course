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
N = 100;
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



vec_d = 0:l/100:3*l-l/100;
c1 = zeros(1,length(vec_d));
for idx = 1:length(vec_d)
    d2 = vec_d(idx);
    H = channel_uniform_angle_K(Nt, Nr, t, fd,d1,d2,N,l,K);
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
plot((1:length(vec_d))/100,c1)
xlabel('spacing (wavelengths)')
ylabel('correlation')
title("Spatial correlation vs antenna spacing")




%% 1 a).

H = channel_uniform_angle_K(Nt, Nr, t, fd,d1,d2,N,l,K);
rx1 = H(:,1,1);
rx2 = H(:,2,1);
figure;
semilogy(t,abs(rx1))
hold on
semilogy(t,abs(rx2))
title("1(a)")
xlabel('time (s)')
ylabel('signal amplitude (db)')
title("Signal amplitude on two antennas")
legend("rx1","rx2")
ylim([1e-4,1e1])
xlim([0, T])

%% 1 b).
T = 10;
%time axis
t = 0:1/fs:T-1/fs;
Niter = 100;
P1 = zeros(1,Niter);
P2 = zeros(1,Niter);
for i = 1:Niter
H = channel_uniform_angle_K(Nt, Nr, t, fd,d1,d2,N,l,K);
rx1 = H(:,1,1);
rx2 = H(:,2,1);
P1(1,i) = sum(abs(rx1).^2)/length(t);
P2(1,i) = sum(abs(rx2).^2)/length(t);
end
figure;
plot(P1)
hold on
plot(P2)
xlabel("Experiment number")
ylabel("Receiver power (W)")
legend("rx1","rx2")
%% 1 c).
T = 0.1;
%time axis
t = 0:1/fs:T-1/fs;

temp_ac = zeros(1,length(t));
parfor tau = 1:length(t)
    for n = 1:10
        H = channel_uniform_angle_K(Nt, Nr, t, fd,d1,d2,N,l,K);
        rx1 = H(:,1,1);
        temp_ac(tau) = temp_ac(tau) + rx1(1)*rx1(tau)';
    end
    temp_ac(tau) = real(temp_ac(tau)/100);
end
figure;
plot(t,temp_ac/max(temp_ac),'-o')
hold on
plot(t(1:100),besselj(0,2*pi*fd*t(1:100)),'X')
xlabel("time (s)")
ylabel("Correlation")
title("Temporal correlation across 1 antenna fd = 100hz")
legend("simulated","theoretical distribution")



%% 1 d).
%Time duration
T = 0.1;
%time axis
t = 0:1/fs:T-1/fs;
%max doppler
fd = 15;

H = channel_uniform_angle_K(Nt, Nr, t, fd,d1,d2,N,l,K);
rx1 = H(:,1,1);
rx2 = H(:,2,1);
figure;
semilogy(t,abs(rx1))
hold on
semilogy(t,abs(rx2))
title("1(a)")
xlabel('time (s)')
ylabel('signal amplitude (db)')
title("Signal amplitude on two antennas fd=15 hz")
legend("rx1","rx2")
ylim([1e-4,1e1])
xlim([0, T])


T = 10;
%time axis
t = 0:1/fs:T-1/fs;
H = channel_uniform_angle_K(Nt, Nr, t, fd,d1,d2,N,l,K);
rx1 = H(:,1,1);
rx2 = H(:,2,1);
sum(abs(rx1).^2)/length(t)
sum(abs(rx2).^2)/length(t)


T = 0.1;
%time axis
t = 0:1/fs:T-1/fs;

temp_ac = zeros(1,length(t));
parfor tau = 1:length(t)
    for n = 1:10
        H = channel_uniform_angle_K(Nt, Nr, t, fd,d1,d2,N,l,K);
        rx1 = H(:,1,1);
        temp_ac(tau) = temp_ac(tau) + rx1(1)*rx1(tau)';
    end
    temp_ac(tau) = real(temp_ac(tau)/100);
end
figure
plot(t,temp_ac/max(temp_ac),'-o')
hold on
xlabel("time (s)")
ylabel("Correlation")
title("Temporal correlation across 1 antenna fd = 15 hz")

plot(t(1:100),besselj(0,2*pi*fd*t(1:100)),'X')
legend("simulated","theoretical distribution")


