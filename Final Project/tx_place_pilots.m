function [tx1_complex_sym,tx2_complex_sym, tx_pilots_tx1, tx_pilots_tx2] = tx_place_pilots(tx1_complex_sym,tx2_complex_sym, pilot_scs_tx1_idx, pilot_scs_tx2_idx)
Np = length(pilot_scs_tx2_idx);
%assuming QPSK pilots
pilot_bits = floor(2*rand(2*Np,1));
pilots_dec = bit2int(reshape(pilot_bits,[2,Np]),2);
tx_pilots_tx1 = qammod(pilots_dec,2^2,'UnitAveragePower',true);
tx1_complex_sym(pilot_scs_tx1_idx) = tx_pilots_tx1;
%null tx2 pilot scs on tx1 
tx1_complex_sym(pilot_scs_tx2_idx) = 0;


pilot_bits = floor(2*rand(2*Np,1));
pilots_dec = bit2int(reshape(pilot_bits,[2,Np]),2);
tx_pilots_tx2 = qammod(pilots_dec,2^2,'UnitAveragePower',true);
tx2_complex_sym(pilot_scs_tx2_idx) = tx_pilots_tx2;
%null tx2 pilot scs on tx1
tx2_complex_sym(pilot_scs_tx1_idx) = 0;

end