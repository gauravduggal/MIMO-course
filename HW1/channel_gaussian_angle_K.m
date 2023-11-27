
function [H] = channel_gaussian_angle_K(Nt, Nr, t, fd, d1,d2,N,l,K,u,sigma)
H = zeros(length(t),Nr,Nt);
phinvec = 2*pi*rand(1,N-1);
thetavec = normrnd(u*pi/180,sigma*pi/180,[1,N-1]);
omegavec = normrnd(u*pi/180,sigma*pi/180,[1,N-1]);
phin0 = 2*pi*rand(1);  
    for idrx = 1:Nr
        for idtx = 1:Nt
%                 phin0 = 2*pi*rand(1);
                
                %dominant path at angle = 0
                thetan = 0;
                omegan = 0;
%                 p_dominant = zeros(1,length(t));
                p_dominant = exp(1j*phin0)*exp(1j*2*pi*fd*cos(thetan)*t)*exp(-1j*2*pi/l*d1*(idtx-1)*sin(omegan))*exp(-1j*2*pi/l*d2*(idrx-1)*sin(thetan));
                p_dominant = p_dominant/(std(p_dominant));
                mpc_combined = zeros(1,length(t));
                for idn = 1:N-1
                    phin = phinvec(1,idn);
                    %AoA
                    thetan = thetavec(1,idn);
                    %AoD
                    omegan = omegavec(1,idn);
                    temporal = exp(1j*2*pi*fd*cos(thetan)*t+phin);
%                     mpc_phase = exp(1j*phin);
                    tx_side = exp(-1j*2*pi/l*d1*(idtx-1)*sin(omegan));
                    rx_side = exp(-1j*2*pi/l*d2*(idrx-1)*sin(thetan));
                    mpc = temporal*tx_side*rx_side;
                    
                    mpc_combined = mpc_combined + sqrt(1/(N-1))*mpc; 
                end
                 mpc_combined = mpc_combined/(std(mpc_combined));
                H(:,idrx, idtx) = sqrt((K)/(K+1))*transpose(p_dominant(1,:)) + sqrt((1)/(K+1))*transpose(mpc_combined);
%               
        end
            
    end
%     for idrx = 1:Nr
%         for idtx = 1:Nt
%             rms_amp = std(H(:,idrx,idtx));
%             H(:,idrx,idtx) = H(:,idrx,idtx)/rms_amp;
%         end
%     end
    
end
