function [m,v,ph]=Division_Alpha_Gamma_D(Mnl , Snl ,Xci_p0,m_alpha,v_alpha,Poid,D)
pdf_Gauss_log=@(t,m,s2,d) (-0.5*(t-m )'*inv(s2)*(t-m ));
% ss=diag(Snl);
    mnl=[[(Mnl(:,1)-2*pi) (Mnl(:,2)-2*pi)];[(Mnl(:,1)-2*pi) Mnl(:,2)];[Mnl(:,1) (Mnl(:,2)-2*pi)];  Mnl  ;[(Mnl(:,1)+2*pi) Mnl(:,2)];[Mnl(:,1) (Mnl(:,2)+2*pi)]; [(Mnl(:,1)+2*pi) (Mnl(:,2)+2*pi)]; [(Mnl(:,1)+2*pi) (Mnl(:,2)-2*pi)]; [(Mnl(:,1)-2*pi) (Mnl(:,2)+2*pi)]];
    snl=[Snl;Snl;Snl;Snl;Snl;Snl;Snl;Snl;Snl];
    Xci=[Xci_p0;Xci_p0;Xci_p0;Xci_p0;Xci_p0;Xci_p0;Xci_p0;Xci_p0;Xci_p0];
MM=size(mnl,1);
ph=[];
for i=1:MM    
ph=[ph;(Xci(i)+0.5*log(det(snl((i-1)*D+1:i*D,:))./det(v_alpha))-(pdf_Gauss_log( inv(v_alpha)*m_alpha.',inv(snl((i-1)*D+1:i*D,:))*mnl(i,:).',inv(v_alpha)-inv(snl((i-1)*D+1:i*D,:)),D ))' - (pdf_Gauss_log( mnl(i,:).',0,v_alpha ,D ))' + (pdf_Gauss_log( m_alpha.',0,snl((i-1)*D+1:i*D,:),D ))')];
end
[val,ind]=max(ph);
q=[(ind-1)*D+1 ind*D]';
q=q(:);
Smax=snl(q,:);
Mmax=mnl(ind,:);
ph=ph(ind);        
% S =v_alpha - Smax;
% m =(Smax*inv(S)*m_alpha.' - v_alpha*inv(S)*Mmax.').';
v = inv((inv(v_alpha) - inv(Smax))) ;
m=(v*(inv(v_alpha)*m_alpha.' - inv(Smax)*Mmax.')).';
while ( (m(1)> pi) || (m(1)< -pi) || (m(2)> pi) || (m(2)< -pi) )
    m(m<-pi)=m(m<-pi)+2*pi;
    m(m>pi)=m(m>pi)-2*pi;
end
end