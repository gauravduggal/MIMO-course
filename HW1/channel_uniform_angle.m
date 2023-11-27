
function [H] = channel_uniform_angle(Nt, Nr, t, fd, d1,d2,N,l)

H = zeros(Nr,Nt,length(t));
    for idrx=1:Nr
        for idtx = 1:Nt
            for n=1:N
                phin = 2*pi*rand(1);
                %AoA
                thetan = 2*pi*rand(1);
                %AoD
                omegan = 2*pi*rand(1);
                %NLOS path at uniform spread angles
                tx_side = exp(-1j*2*pi/l*d1*(idtx-1)*sin(omegan));
                rx_side = exp(-1j*2*pi/l*d2*(idrx-1)*sin(thetan));
                for idt = 1:length(t)
                    temporal = exp(1j*2*pi*fd*cos(thetan)*t(idt))*exp(1j*2*pi*phin);
                    H(idrx,idtx,idt) = H(idrx,idtx,idt) + temporal*tx_side*rx_side;
                end

            end
        end
    end
end