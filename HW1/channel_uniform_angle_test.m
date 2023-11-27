
function [H] = channel_uniform_angle_test(Nt, Nr, t, fd, d1,d2,N,l)

H = zeros(Nr,Nt,length(t));
    for idrx=1:Nr
        for idtx = 1:Nt
            for idt = 1:length(t)
                for n=1:N
                    %AoA and AoD at uniform spread angles
                    phin = 2*pi*rand(1);
                    %AoA
                    thetan = 2*pi*rand(1);
                    %AoD
                    omegan = 2*pi*rand(1);
                    tx_side = exp(-1j*2*pi/l*d1*(idtx-1)*sin(omegan));
                    rx_side = exp(-1j*2*pi/l*d2*(idrx-1)*sin(thetan));
                    %random phase based on path length
                    temporal = exp(1j*2*pi*fd*cos(thetan)*t(idt))*exp(1j*phin);
                    mpc = temporal*tx_side*rx_side;
                    H(idrx,idtx,idt) = H(idrx,idtx,idt) + mpc;
                end
            H(idrx,idtx,idt) = H(idrx,idtx,idt) /abs(H(idrx,idtx,idt));    
            end
            
        end
    end
end