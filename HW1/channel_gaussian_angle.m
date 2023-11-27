function [H] = channel_gaussian_angle(Nt, Nr, t, fd, u, sigma,d1,d2,N,l)
H = zeros(Nr,Nt,length(t));
    for i=1:Nr
        for j = 1:Nt
            %Dominant LOS path at zero degrees
            thetan = 0;
            omegan = 0;
            phin = 2*pi*rand(1);
            for k = 1:length(t)
                H(i,j,k) = 1/sqrt(2)*exp(1j*2*pi*fd*cos(thetan)*t(k)+phin)*exp(-1j*2*pi/l*d1*(j-1)*sin(omegan))*exp(-1j*2*pi/l*d2*(i-1)*sin(thetan));
                %NLOS path at gaussian spread angles
                for n = 1:N-1
                    %phase due to path propagation
                    phin = 2*pi*rand(1);
                    %AoA
                    thetan = normrnd( u*pi/180 , sigma*pi/180);%2*pi*rand(1);
                    %AoD
                    omegan = normrnd( u*pi/180 , sigma*pi/180);
                    H(i,j,k) = H(i,j,k) + 1/sqrt(2*(N-1))*exp(1j*2*pi*fd*cos(thetan)*t(k)+phin)*exp(-1j*2*pi/l*d1*(j-1)*sin(omegan))*exp(-1j*2*pi/l*d2*(i-1)*sin(thetan));
                end
            end
            
        end
    end
    return
end

