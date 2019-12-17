function [cout,PER]=Decodeur(hDec,Data,LLR,Nc)
Nb=size(Data,1);
LLR = LLR';
LLR = LLR(:);
LLR = -LLR(1:Nc);
dec_llr_ph=step(hDec,LLR);
dec_llr_ph=dec_llr_ph(1:Nb);
Binest_ph =  (dec_llr_ph<0) ;
difference     = ( Binest_ph ~= Data);
count =(sum(difference));
PER   = any( difference );

end
