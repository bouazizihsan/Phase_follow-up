function [m,v,ph]=Division_Alpha_Gamma(Mnl , Snl ,Xci_log,m_alpha,v_alpha,poid_log)
pdf_Gauss_log   = @(t,m,s2) ( -abs(t-m ).^2./(2*s2) )- 0.5*log(2*pi*s2);

ph= poid_log - (Xci_log) + pdf_Gauss_log( m_alpha ,Mnl,v_alpha + Snl);
v =1./(1./v_alpha - 1./Snl) ;%(v_alpha .* Snl ./ (v_alpha + Snl));
m = v.*(m_alpha./v_alpha  - Mnl./Snl);% m_alpha./(v_alpha./Snl + 1) + Mnl./(1 + Snl./v_alpha) ;  

if(v<0)
 stp=1;   
end
end