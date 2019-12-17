function [m,v,ph]=Produit_Alpha_Gamma_D(Mnl , Snl ,Xci_p0,m_alpha,v_alpha,Poid,D)
%pdf_Gauss_log=@(t,m,s2,d) ((-0.5*(Operator(t-m ))'*inv(s2)*(Operator(t-m ))) - log( (2*pi).^(d/2)*det(s2).^(0.5) ));
pdf_Gauss_log=@(t,m,s2,d) ((-0.5*((t-m ))'*inv(s2)*((t-m ))) - log( (2*pi).^(d/2)*det(s2).^(0.5) ));

ss=diag(Snl);
if(ss==Inf)
    ph=Xci_p0 + Poid;
    m= m_alpha;
    v= v_alpha;
else
%     mnl=Mnl; % [[(Mnl(:,1)-2*pi) (Mnl(:,2)-2*pi)];[(Mnl(:,1)-2*pi) Mnl(:,2)];[Mnl(:,1) (Mnl(:,2)-2*pi)];  Mnl  ;[(Mnl(:,1)+2*pi) Mnl(:,2)];[Mnl(:,1) (Mnl(:,2)+2*pi)]; [(Mnl(:,1)+2*pi) (Mnl(:,2)+2*pi)]; [(Mnl(:,1)+2*pi) (Mnl(:,2)-2*pi)]; [(Mnl(:,1)-2*pi) (Mnl(:,2)+2*pi)]]; %
%     snl=Snl; % [Snl;Snl;Snl;Snl;Snl;Snl;Snl;Snl;Snl];%
%     Xci=Xci_p0; % [Xci_p0;Xci_p0;Xci_p0;Xci_p0;Xci_p0;Xci_p0;Xci_p0;Xci_p0;Xci_p0];%
%     MM=size(mnl,1);
%     ph=[];
%     for i=1:MM
%         ph=[ph;Xci(i) + Poid + (pdf_Gauss_log( m_alpha.' ,mnl(i,:).',v_alpha + snl((i-1)*D+1:i*D,:),D ))'];
%     end
%     [val,ind]=max(ph);
%     q=[(ind-1)*D+1 ind*D]';
%     q=q(:);
%     Smax=snl(q,:);
%     Mmax=mnl(ind,:);
%     ph=ph(ind);
%     S =v_alpha + Smax;

    ph=Xci_p0 + Poid + (pdf_Gauss_log( m_alpha.' ,Mnl.',v_alpha + Snl,D ))';
    v = inv(inv(v_alpha) + inv(Snl)); % inv(inv(v_alpha) + inv(Smax)) ;
    %m = Operator(m_alpha+(v*inv(Snl)*(Operator(Mnl-m_alpha)).').'); %(Smax*inv(S)*m_alpha.' + v_alpha*inv(S)*Mmax.').';
    m = (m_alpha+(v*inv(Snl)*((Mnl-m_alpha)).').'); 
end

end