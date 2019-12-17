function [m,v,ph]=Produit_Alpha_Gamma(Mnl , Snl ,Xci_log,m_alpha,v_alpha,poid_log)

pdf_Gauss_log   = @(t,m,s2) ( -abs(Operator(t-m )).^2./(2*s2) )- 0.5*log(2*pi*s2);
% pdf_Gauss_log   = @(t,m,s2) ( -abs((t-m )).^2./(2*s2) )- 0.5*log(2*pi*s2);

if((Snl==Inf))
    ph= (Xci_log) + poid_log;
else
    ph= (Xci_log)+ poid_log + pdf_Gauss_log( m_alpha ,Mnl,v_alpha + Snl);
end
v =1./(1./v_alpha + 1./Snl) ;%(v_alpha .* Snl ./ (v_alpha + Snl));
% m = ((m_alpha.*Snl + Mnl.*v_alpha)./(v_alpha + Snl));% m_alpha./(v_alpha./Snl + 1) + Mnl./(1 + Snl./v_alpha) ;
% m= v.*(m_alpha./v_alpha + Mnl./Snl );
 m= (m_alpha + v.*(Operator(Mnl-m_alpha))./Snl) ;

end