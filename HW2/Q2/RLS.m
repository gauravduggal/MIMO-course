close all
clear cll
clc

Ns = 300;

M = 1;
N = 4;
c = 3e8;
f = 9e9;
lambda = c/f;
d = lambda/2;



AoA = 45;
AoD = 0;



%AoA for interferer
AoA_int = 25;
AoA_int2 = 0;
a2 = array_response(AoA_int2,N,lambda,d);
a1 = array_response(AoA_int,N,lambda,d);


chi = 0.99;
K = zeros(N,Ns);
e = zeros(1,Ns);
w_rls = zeros(N,Ns);
w_mmse = zeros(N,Ns);
P = zeros(N,N,Ns);
P(:,:,1) = eye(N);

rx_frame = zeros(N,Ns);
tx_frame = zeros(Ns,1);
data = zeros(Ns,1);
iin = zeros(Ns,1);
iin2 = zeros(Ns,1);
for idx = 1:Ns
    data(idx,1) = randi([0,3],[1,1]);
    iin(idx,1) = randi([0,3],[1,1]);
    iin2(idx,1) = randi([0,3],[1,1]);
    iin_sym =  1/sqrt(2)*qammod(iin(idx,1),4);
    iin_sym2 =  1/sqrt(2)*qammod(iin2(idx,1),4);


    % modulate the data symbols
    syms = 1/sqrt(2)*qammod(data(idx,1),4);


    %     %transmit symbols along AoD
    tx_frame(idx,1) = transmit(syms,ones(M,1),AoD,lambda,d);


    AoA_soi = AoA + 10*sin(2*pi*Ns/10*idx/Ns);
    a0 = array_response(AoA_soi,N,lambda,d);
    % NxNs receive matrix
    rx_frame(:,idx) = a0*transpose(tx_frame(idx,1))+1/sqrt(2)*(a1*transpose(iin_sym)+a2*transpose(iin_sym2));

    snr_db = 10;% SNR_vec(idy);


    %add noise
    [rx_frame(:,idx),ns,txx] = add_noise(rx_frame(:,idx),snr_db);
end


for k = 2:Ns


    % RLS
    K(:,k) = (chi^-1*P(:,:,k-1)*rx_frame(:,k))/(1+chi^-1*rx_frame(:,k)'*P(:,:,k-1)*rx_frame(:,k));
    e(k) = tx_frame(k,1) - w_rls(:,k-1)'*rx_frame(:,k);
    w_rls(:,k) = w_rls(:,k-1)+K(:,k)*e(k)';
    P(:,k) = chi^-1*P(:,k-1) - chi^-1*K(:,k)*rx_frame(:,k)'*P(:,k-1);

    % MMSE
    Rrr = rx_frame(:,1:k)*rx_frame(:,1:k)';
    w_mmse(:,k) = pinv(Rrr)*rx_frame(:,1:k)*conj(tx_frame(1:k,1))/k;

end


% w_mmse = (Rrr)^-1*rx_pilots*tx_pilots';

beampattern(w(:,2),N,lambda,d,true);
beampattern(w(:,3),N,lambda,d,true);
beampattern(w(:,4),N,lambda,d,true);
beampattern(w(:,10),N,lambda,d,true);
beampattern(w(:,Ns),N,lambda,d,true);
beampattern(w_mmse(:,Ns),N,lambda,d,true);



Ntrials = 1000;
ber_mmse = zeros(1,Ns);
ber_rls = zeros(1,Ns);

parfor k = 2:Ns
    for n = 1:Ntrials
        %generate Npilots and (Ns-Npilots) data symbols
        data = randi([0,3],[Ns,1]);
        iin = randi([0,3],[Ns,1]);
        iin2 = randi([0,3],[Ns,1]);
        iin_sym =  1/sqrt(2)*qammod(iin,4);
        iin_sym2 =  1/sqrt(2)*qammod(iin2,4);


        % modulate the data symbols
        syms = 1/sqrt(2)*qammod(data,4);


        %     %transmit symbols along AoD
        %     w = ones(M,1);
        tx_frame = transmit(syms,ones(M,1),AoD,lambda,d);


        AoA_soi = AoA + 10*sin(2*pi*Ns/10*idx/Ns);
        a0 = array_response(AoA_soi,N,lambda,d);


        % NxNs receive matrix
        rx_frame = a0*transpose(tx_frame)+1/sqrt(2)*(a1*transpose(iin_sym)+a2*transpose(iin_sym2));
        snr_db = 5;% SNR_vec(idy);


        %add noise
        [rx_frame,ns,txx] = add_noise(rx_frame,snr_db);



        z_data_RLS = w(:,k)'*rx_frame(:,1:Ns);
        z_data_MMSE = w_mmse(:,k)'*rx_frame(:,1:Ns);

        data_rx_rls = qamdemod(z_data_RLS,4)';
        data_rx_mmse = qamdemod(z_data_MMSE,4)';
        ber_mmse(k) = ber_mmse(k) + 1/Ntrials*(100-100*sum(data==data_rx_mmse)/(Ns));
        ber_rls(k) = ber_rls(k) + 1/Ntrials*(100-100*sum(data==data_rx_rls)/(Ns));
    end
end

plot(2:150,ber_mmse(2:150));
hold on
plot(2:150,ber_rls(2:150));
xlabel("time axis (symbols)")
ylabel("Bit Error Rate")
title("RLS vs MMSE")
grid on
legend("MMSE","RLS")

figure;
AoA_soi = AoA + 10*sin(2*pi*Ns/10*(1:Ns)/Ns);
plot(1:100,AoA_soi(1:100))
xlabel("time (symbols)")
ylabel("AoA of Signal of interest")
grid on