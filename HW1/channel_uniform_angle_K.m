
function [H] = channel_uniform_angle_K(Nt, Nr, t, fd, d1,d2,N,l,K)
H = zeros(Nr,Nt,length(t));
    for idrx = 1:Nr
        for idtx = 1:Nt
                phin0 = 2*pi*rand(1);
                %dominant path at angle = 0
                thetan = 0;
                omegan = 0;
                p_dominant = zeros(1,length(t));
                for idt = 1:length(t)
                    p_dominant(idt) = exp(1j*2*pi*fd*cos(thetan)*t(idt)+phin0)*exp(-1j*2*pi/l*d1*(idtx-1)*sin(omegan))*exp(-1j*2*pi/l*d2*(idrx-1)*sin(thetan));
                end
                temp = zeros(1,length(t));
                for idn = 1:N-1
                    phin = 2*pi*rand(1);
                    %AoA
                    thetan = 2*pi*rand(1);
                    %AoD
                    omegan = 2*pi*rand(1);
                    
                    for idt = 1:length(t)
                        temp = temp + exp(1j*2*pi*fd*cos(thetan)*t(idt)+phin)*exp(-1j*2*pi/l*d1*(idtx-1)*sin(omegan))*exp(-1j*2*pi/l*d2*(idrx-1)*sin(thetan));
                    end
                end
                H(idrx,idtx,idt) = sqrt(K/(K+1))*p_dominant + sqrt(1/(K+1))*1/(N-1)*temp; 
            end
            
        end
end
