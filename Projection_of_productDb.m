function [Mnl,S]=Projection_of_productDb(Mnl_p0,Snl,Xci_p0,parameters1,parameters2,sigw2,D,M1,M2,Phase_noise)

pdf_Gauss=@(t,m,s2,d) exp(-0.5*sum(((t-m )'*inv(s2).*(t-m ).'),2) )./((2*pi).^(d/2)*det(s2).^(0.5));
M11=M1(:);
M22=M2(:);
M=size(Mnl_p0,1);

PHI_nl_len_p0=[];
Mnl_tilde_p0=[];
Snl_tilde=[];

if(diag(parameters2)==Inf)%(diag(alpha_p0{n-1,2})==(ones(D,1)*Inf))
    MM=1; % need to plot fig
    PHI_nl_len_p0=Xci_p0-log((2*pi));
    Mnl_tilde_p0 = Mnl_p0;
    Snl_tilde    =  Snl;
else
    mnl=[[(Mnl_p0(:,1)-2*pi) (Mnl_p0(:,2)-2*pi)];[(Mnl_p0(:,1)-2*pi) Mnl_p0(:,2)];[Mnl_p0(:,1) (Mnl_p0(:,2)-2*pi)];  Mnl_p0  ;[(Mnl_p0(:,1)+2*pi) Mnl_p0(:,2)];[Mnl_p0(:,1) (Mnl_p0(:,2)+2*pi)]; [(Mnl_p0(:,1)+2*pi) (Mnl_p0(:,2)+2*pi)]; [(Mnl_p0(:,1)+2*pi) (Mnl_p0(:,2)-2*pi)]; [(Mnl_p0(:,1)-2*pi) (Mnl_p0(:,2)+2*pi)]];
    snl=[Snl;Snl;Snl;Snl;Snl;Snl;Snl;Snl;Snl];
    Xci=[Xci_p0;Xci_p0;Xci_p0;Xci_p0;Xci_p0;Xci_p0;Xci_p0;Xci_p0;Xci_p0];
    MM=size(mnl,1);
    for i=1:MM
        [m,v,ph]=Produit_Alpha_Gamma_D(mnl(i,:) , snl((i-1)*D+1:i*D,:) ,Xci(i),parameters1,parameters2,D);
        PHI_nl_len_p0=[PHI_nl_len_p0;ph];
        Mnl_tilde_p0=[Mnl_tilde_p0;m];
        Snl_tilde=[Snl_tilde;v];
    end
    
end
PHI_nl_len_p0=exp(PHI_nl_len_p0-jac_logH(PHI_nl_len_p0.'));%% marginalisation

[val ind] = sort(PHI_nl_len_p0,'descend');

Mnl_tilde_p0=Mnl_tilde_p0(ind(1:M),:);
q=[(ind(1:M)-1)*D+1 ind(1:M)*D]';
q=q(:);

Snl_tilde=Snl_tilde(q,:);
PHI_nl_len_p0=PHI_nl_len_p0(ind(1:M));
PHI_nl_len_p0=PHI_nl_len_p0./sum(PHI_nl_len_p0);


% ind=find(PHI_nl_len_p0<0.001);
% if(ind)
%     Mnl_tilde_p0(ind,:)=[];
%     Snl_tilde(sort([(ind-1)*D+1;ind*D]),:)=[];
%     PHI_nl_len_p0(ind)=[];
%     
% end
% PHI_nl_len_p0=PHI_nl_len_p0./sum(PHI_nl_len_p0);


MM=size(Mnl_tilde_p0,1);

for i=1:MM
    fig=PHI_nl_len_p0(i,:).*pdf_Gauss([M11.';M22.'] ,repmat(Mnl_tilde_p0(i,:).',1,length(M11)) , Snl_tilde((i-1)*D+1:i*D,:),D);
    fig=reshape(fig,[],length(M1));
    figure(5)
    contour(M1,M2,fig);
    hold on

    figure(6)

    mesh(M1,M2,fig);
    hold on
end



Mnl_tilde_p0(Mnl_tilde_p0<-pi)=Mnl_tilde_p0(Mnl_tilde_p0<-pi)+2*pi;
Mnl_tilde_p0(Mnl_tilde_p0>pi)=Mnl_tilde_p0(Mnl_tilde_p0>pi)-2*pi;


[Mnl,S]=GaussianD1(Mnl_tilde_p0,Snl_tilde,PHI_nl_len_p0);% Somme_GaussianD(Mnl_tilde_p0,Snl_tilde,PHI_nl_len_p0) ;

fig=(pdf_Gauss([M11.';M22.'] ,repmat(Mnl.',1,length(M11)),S,D ));
fig=reshape(fig,[],length(M1));
figure(5)
hold on
contour(M1,M2,fig);
plot(Phase_noise(1),Phase_noise(2),'r*');
hold off

figure(6)
mesh(M1,M2,fig);
hold off
S=S+sigw2;
end