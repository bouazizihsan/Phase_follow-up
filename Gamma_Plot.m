function Gamma_Plot(X_Hard,Xci_p0,Mnl_p0,Snl,CC)
Np           = 2^9;
Phase_limit  =pi;
dpp          = 2.*Phase_limit./Np;
p           = linspace(-Phase_limit,Phase_limit,Np)';

pdf_Gauss_log   = @(t,m,s2) ( -abs(t-m ).^2./(2*s2) )-0.5*log((2*pi*s2));
pdf_Gauss   = @(t,m,s2) exp( -abs(t-m ).^2./(2*s2) )./sqrt(2*pi*s2);

gamma=(repmat(Xci_p0.',Np,1).*pdf_Gauss(repmat(p,1,length(X_Hard)),repmat(Mnl_p0.',Np,1),repmat(Snl.',Np,1)));
ggg=gamma(:);
gamma=exp(gamma-jac_logH(ggg.'));
Sgamma=sum(gamma,2);
Sgamma=Sgamma./sum(Sgamma);

figure(1)
plot(p,Sgamma(:,1),CC);

end