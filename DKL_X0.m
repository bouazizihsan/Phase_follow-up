function [Result,Coefs] = DKL_X0(Moy,S,PH,NClass,Threshold,Method)
% clc;
% NClass=2;
% epsilon=2;
Np=2^9;
pdf_Gauss   = @(t,m,s2) exp( -abs(t-m ).^2./(2*s2) )./sqrt(2*pi*s2);
Phase_limit = pi;
pp          = linspace(-Phase_limit,Phase_limit,Np)';

% % PH=[0.5;0.15];
% %
% % S=[1.5;0.35];
% % for i=1:10
% % Moy=[0; pi/6+pi/10*(i-1)];

% % graph=repmat(PH.',Np,1).*pdf_Gauss(repmat(pp,1,length(PH)),repmat(Moy.',Np,1),repmat(S.',Np,1));
% % figure(1)
% % plot(pp,graph)

%___Methode de chercher le point d'intersection des deux gaussiannes ______

phi=PH;
m=Moy;
sig=S;

Result=[];
Coefs=[];
j=1;
while ((j <= NClass-1) && (norm(phi)>0))
    [vph,ind]=max(phi);
    phi(ind)=[];
    vg=sig(ind);
    sig(ind)=[];
    vm=m(ind);
    m(ind)=[];
    
    [x0,indx]= ChercheX0(vg,sig,vm,m,vph./sqrt(2*pi*vg),phi./sqrt(2*pi*sig));
    
    ind = find( [(vph.*pdf_Gauss(x0,vm,vg) + phi(indx).*pdf_Gauss(x0,m(indx),sig(indx))) > (phi(indx).*pdf_Gauss(m(indx),m(indx),sig(indx))+ vph.*pdf_Gauss(m(indx),vm,vg))] + [(vph.*pdf_Gauss(x0,vm,vg)+ phi(indx).*pdf_Gauss(x0,m(indx),sig(indx))) > (phi(indx).*pdf_Gauss(vm,m(indx),sig(indx))+ vph.*pdf_Gauss(vm,vm,vg))] + [ vph.*pdf_Gauss(vm,vm,vg) > phi(indx).*pdf_Gauss(vm,m(indx),sig(indx)) ] );
    fact=[vph;phi(ind)];
    u=[vm;m(ind)];
    var=[vg;sig(ind)];
    switch Method
        case 'Gauss'
            
            [mt,st]=Gaussian1(u,var,fact);
            if(st<Threshold)
                Result= [Result [mt ;st;sum(fact)]];
            else
                Result= [Result [mt ;Inf;sum(fact)]];
                
            end
        case 'Tikhonov'
            
            [Z]=CMVM(fact,(1./var).*exp(1j*u),Threshold);
            Coefs=[Coefs;sum(fact)];
            Result= [Result;Z];
            mt=angle(Z);
            st=1./abs(Z);
            
    end
       gauss= sum(fact)*pdf_Gauss(pp, mt,st);
       
    phi(ind)=[];
    sig(ind)=[];
    m(ind)  =[];
    j=j+1;
end

% % figure(2)
% % plot(pp,gauss,'b')
% % hold on
% % if(norm(phi))
% % plot(pp,(repmat(phi.',Np,1).*pdf_Gauss(repmat(pp,1,length(phi.')),repmat(m.',Np,1),repmat(sig.',Np,1))),'--')
% % end  
hold off
if(norm(phi))
    switch Method
        case 'Gauss'
            [ mt,st]=Gaussian1(m,sig,phi);
            
            if(st<Threshold)
                Result= [Result [mt ;st;sum(phi)]];
                
            else
                Result= [Result [mt ;Inf;sum(phi)]];
            end
        case 'Tikhonov'
            
            [Z]=CMVM(phi,(1./sig).*exp(1j*m),Threshold);
            Coefs=[Coefs;sum(phi)];
            Result= [Result;Z];
            
            
    end
        
end

%___ Methode de Kullback Leiber Divergence ______
% % phi=PH;
% % m=Moy;
% % sig=S;
% % j=1;
% % gauss2=[];
% % while ((j <= NClass-1) && (norm(phi)>0))
% %     [vph,indx]=max(phi);
% %     f=DKL(m(indx),m,sig(indx),sig);
% %     f
% %     aux=find(f<epsilon);
% %     sum_phi=sum(phi(aux));
% %     [MT,ST]=Gaussian1(m(aux).',sig(aux).',phi(aux).');
% %     gauss2=[gauss2  (sum_phi*pdf_Gauss(pp,MT,ST))];
% %     phi(aux)=[];
% %     sig(aux)=[];
% %     m(aux)  =[];
% %     j=j+1;
% % end
% %
% % hold
% % hold on
% % plot(pp,gauss2,'g')
% % if(norm(phi))
% % plot(pp,repmat(phi,Np,1).*pdf_Gauss(repmat(pp,1,length(phi)),repmat(m,Np,1),repmat(sig,Np,1)),'k--')
% % end
% % plot(pp,(phi*pdf_Gauss(repmat(pp,1,length(phi)),repmat(m,Np,1),repmat(sig,Np,1))')' + gauss2,'k->')
% % plot(pp,(PH*pdf_Gauss(repmat(pp,1,length(PH)),repmat(Moy,Np,1),repmat(S,Np,1)).').', 'c--' )
% % legend('Appx1','reste1','Appx1+reste1','AppxKL','reste2','AppxKL+reste2','somme')
% % hold off

end