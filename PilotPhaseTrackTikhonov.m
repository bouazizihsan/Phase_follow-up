function [BCJR]=PilotPhaseTrackTikhonov( y, N0, sigw2, X, Lframe)
% pdf_Tikh=@(t,k,u) exp(k.*cos(t-u))./(2*pi.*besseli(0,k));
pdf_Tikh=@(t,z) exp(real(z.*exp(-1j*t)))./(2*pi.*besseli(0,abs(z)));
pdf_Gauss   = @(t,m,s2) exp( -abs(t-m ).^2./(2*s2) )./sqrt(2*pi*s2);

Np           = 2^9;
pp           = linspace(-pi,pi,Np)'; 

D=1;
Ns           =length(y);
Z_alpha=zeros(Ns,D);

for n=2:Ns
    
    X_Hard=X(n-1);    
    z=Z_alpha(n-1) + (2/N0)*y(n-1)*conj(X_Hard);
    
% %     g=ifft(fft(lambda.*pdf_Tikh(pp,abs(z),angle(z))).*(fft(pdf_Tikh(pp,1./sigw2,0))));
%     g=ifft( fft(pdf_Tikh(pp,z)).*(fft(fftshift(pdf_Tikh(pp,(1./sigw2).*exp(1j*0))))).^Lframe );
%     g=g./sum(g);
%     z_gamma =gamma_fct(sigw2*Lframe,z); 
%     %     gg=pdf_Tikh(pp,abs(zT),angle(zT));
%     gg=pdf_Tikh(pp,z_gamma);
%     gg=gg./sum(gg);
%     figure;
%     plot(pp,g,'g');
%     hold on
%     plot(pp,gg,'r');
%     hold off
    z_gamma =gamma_fct(sigw2,z); 
    Z_alpha(n)=z_gamma;
    
end

%Beta
Z_beta=zeros(Ns,D);

for n=Ns-1:-1:1
    
    X_Hard=X(n+1);

    z=Z_beta(n+1) + (2/N0)*y(n+1)*conj(X_Hard);
      
%    %     g=ifft(fft(lambda.*pdf_Tikh(pp,abs(z),angle(z))).*(fft(pdf_Tikh(pp,1./sigw2,0))));
%     g=ifft( fft(pdf_Tikh(pp,z)).*fft(pdf_Gauss(pp,0,sigw2)) );
%     g=g./sum(g);
%     z_gamma =gamma_fct(sigw2,z); 
%     %     gg=pdf_Tikh(pp,abs(zT),angle(zT));
%     gg=pdf_Tikh(pp,z_gamma);
%     gg=gg./sum(gg);
%     plot(pp,g,'g');
%     hold on
%     plot(pp,gg,'r');
%     hold off
    z_gamma =gamma_fct(sigw2,z); 
    Z_beta(n)=z_gamma;
       
end

BCJR.Z_beta=Z_beta.'; 
BCJR.Z_alpha=Z_alpha.';
end
function z=gamma_fct(sig,z)

z=z./(1+abs(z).*sig);

end

% function z=gamma_fct(sig,z)
% k1=1/sig;
% k2=abs(z);
% u=angle(z);
% circular_var2=1-(besseli(1,k1)./besseli(0,k1) .* besseli(1,k2)./besseli(0,k2) );
% k=1/(2*circular_var2);
% 
% z=k.*exp(1j*u);
% 
% end