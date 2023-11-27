function [Pb] = BER_Maximal_ratio_transmission_2rx_Ltx(L,c)
%Performance of Maximal Ratio Transmission with Two Receive Antennas


temp1 = 0;
for i=0:L-2
    temp1 = temp1+gcna(c,i+L+1,2)*nchoosek(i+L,L);
end

temp2 = 0;
for i=0:L-1
    temp2 = temp2+gcna(c,i+L,2)*nchoosek(i+L-1,L-1);
end

temp3 = 0;
for i=0:L
    temp3 = temp3+gcna(c,i+L-1,2)*nchoosek(i+L-2,L-2);
end

Pb = L*gcna(c,L+1,1)-2*(L-1)*gcna(c,L,1) ...
+L*gcna(c,L-1,1) ...
-L*temp1 ...
+2*(L-1)*temp2 ...
-L*temp3;
end