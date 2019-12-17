function [Mnl,Sigma,PHI,Pc]=Projection_of_product(Mnl_p0,Snl,Xci_p0,parameters1,parameters2,parameters3,InfP,is_BCRP,etape,sigw2,Phase_noise)
NClass=4;

M=size(Mnl_p0,1);
ph=[];
mm=[];
vv=[];

ph=[];
mm=[];
vv=[];
for j=1:size(parameters3,1)
    if(parameters2(j)==Inf)
        phi= (Xci_p0) + parameters3(j);
        mmm = Mnl_p0;
        vvv = Snl;
    else
        mnl=Mnl_p0;
        snl=Snl;
        Xci=Xci_p0;
        %         mnl=[(Mnl_p0-2*pi) Mnl_p0 (Mnl_p0+2*pi)];
        %         mnl=mnl.';
        %         mnl=mnl(:);
        %         snl=[Snl Snl Snl];
        %         snl=snl.';
        %         snl=snl(:);
        %         Xci=[Xci_p0 Xci_p0 Xci_p0];
        %         Xci=Xci.';
        %         Xci=Xci(:);
        
        [mmm,vvv,phi]=Produit_Alpha_Gamma(mnl , snl ,Xci,parameters1(j),parameters2(j),parameters3(j));
    end
    ph=[ph;phi];
    mm=[mm;mmm];
    vv=[vv;vvv];
end
ph=exp(ph -jac_logH(ph'));
m=mm;
v=vv;
Pc.v=v;
Pc.m=m;
Pc.ph=ph;
% if(length(ph)>2)
%     [ph,ind]=max(reshape(ph,3,[]));
%     ph=ph';
%     mm=reshape(mm,3,[]);
%     vv=reshape(vv,3,[]);
%
%     for i=1:size(mm,2)
%         m(i,1)=mm(ind(i),i);
%         v(i,1)=vv(ind(i),i);
%     end
% else
%     m=mm;
%     v=vv;
% end

% % [val ind] = sort(ph,'descend');
% % m=m(ind(1:M));
% % v=v(ind(1:M));
% % ph=ph(ind(1:M));
% % ph=ph./sum(ph);


% m(m<-pi)=m(m<-pi)+2*pi;
% m(m>pi)=m(m>pi)-2*pi;

switch etape
    case 0
        Mnl=m;
        Sigma=v+sigw2;
        PHI=log(ph);
    case 1
        MM=length(ph);
        if(MM>1)
            f=[];
            indx1=[];
            while(length(ph)>NClass)
                
                if(length(f))
                    indice= [repmat(length(ph),length(ph)-1,1),(1:(length(ph)-1))'];
                else
                    indice=(1:length(ph));
                    indice=nchoosek(indice,2);
                    
                end
                
                Mnl=zeros(length(indice),1);
                S=zeros(length(indice),1);
                for j=1:length(indice)
                    q=indice(j,:);
                    q=q(:);
                    [Mnl(j),S(j)]=Gaussian1(m(indice(j,:)),v(q,:),ph(indice(j,:)));
                    
                    f=[f;sum( ph(indice(j,:))./(sum(ph(indice(j,:)))).*DKL(m(indice(j)).',Mnl(j),v(q),S(j)))];
                end
                if(indx1)
                    indice=[indx1;indice];
                    S=[S1;S];
                    Mnl=[Mnl1;Mnl];
                end
                [val ind]=min(f);
                q=(indice(ind,:)).';
                [X,Y]=find((indice==q(1)) + (indice==q(2)));
                f(X)=[];
                
                coef=sum(ph(q));
                v(q)=[];
                m(q)=[];
                ph(q)=[];
                v=[v;S(ind)];
                m=[m;Mnl(ind)];
                ph=[ph;coef];
                S(X)=[];
                S1=S;
                Mnl(X)=[];
                Mnl1=Mnl;
                if(length(ph)>2)
                    indx1=(1:(length(ph)-1)).';
                    indx1=nchoosek(indx1,2);
                else
                    indx1=[];
                end
            end
            Mnl=m;
            Sigma=v+sigw2;
            PHI=log(ph);
            
        else
            PHI=log(ph);
            Sigma=v+sigw2;
            Mnl=m;
            
        end
        %         %%%%%%%         HB-Track
        %         [ind val] = find(ph>epsilon);
        %
        %         m=m(ind,:);
        %
        %         Snl_tilde=Snl_tilde(q,:);
        %         PHI_nl_len_p0=PHI_nl_len_p0(ind,:);
        %         cor=sum(PHI_nl_len_p0);
        %         PHI=PHI_nl_len_p0./sum(PHI_nl_len_p0);
        %         [Mnl,S]=GaussianD1(Mnl_tilde_p0,Snl_tilde,PHI_nl_len_p0);
        %         PHI=0; %log(1)
        
    case 2
        if(M>1)
            %         [Mnl,Sigma]=Gaussian1(m,v,ph);
            %             [Mnl,Sigma]=Somme_Gaussian(m,v,ph);
            
            if(~is_BCRP)
                % ******** if we use F*alpha and Beta(Pilot) to find the best m= argmin(DKL(m,Mnl,v,Sigma) + DKL(m,InfP(1),v,InfP(2))) *********
                [Mnl,Sigma]=Somme_Gaussian(m,v,ph);
                f=(DKL(m,Mnl,v,Sigma) + DKL(m,InfP(1),v,InfP(2)));
                %             [val,ind]=min(f);
                %             Mnl=m(ind);
                %             Sigma=v(ind);
                [val ind] = sort(f);
                [Mnl,Sigma]=Somme_Gaussian(m(ind(1:NClass)),v(ind(1:NClass)),ph(ind(1:NClass))./sum(ph(ind(1:NClass))) );
                
                
                %                 [val,ind]=max(ph);
                %                 Sigma=v(ind);
                %                 Mnl=m(ind);
                
                %                 % ******** if we use F*alpha*Beta(Pilot) *********
                %                 [val ind] = sort(ph);
                %                 [Mnl,Sigma]=Somme_Gaussian(m(ind(1:2)),v(ind(1:2)),ph(ind(1:2))./sum(ph(ind(1:2))));
                
                %                 % ********END if we use F*alpha*Beta(Pilot) *********
                
                %                 % ******** if we use F*alpha + Beta(Pilot) *********
                %                 [Mnl1,Sigma1]=Somme_Gaussian(m,v,ph);
                %                 [Mnl,Sigma]=Somme_Gaussian([Mnl1;InfP(1)],[Sigma1;InfP(2)],[0.5;0.5]);
                %                 Sigma=Sigma.* sqrt(InfP(2).*Sigma1);
                %                 % ********END if we use F*alpha+Beta(Pilot) *********
            else
                [Mnl,Sigma]=Somme_Gaussian(m,v,ph);
                
            end
        else
            
            Mnl=m;
            Sigma=v;
        end
        %         Sigma=Sigma+sigw2;
        PHI=0;
end
% Taw=sum(TawX.*exp(parameters3)./(2*pi) + TawX.*TawAlpha + TawAlpha.*exp(Xci_p0)./(2*pi)) ;

end
function f=DKL(A,B,C,D)
%f= sum(A.*log(A./repmat(B,1,size(A,2))),1);
dist=abs(B-A);
[ind v]=find(dist>pi);
dist(ind)=dist(ind)-2*pi;

f=dist.^2./(2*C) + D./(2*C) +0.5*log(C./D) - 0.5;

end
% function f=DKL(u1,u2,C,D)
% %f= sum(A.*log(A./repmat(B,1,size(A,2))),1);
% dist=abs(u1-u2);
% [ind v]=find(dist>pi);
% dist(ind)=dist(ind)-2*pi;
%
% f=0.5*log(D./C) + (C+dist.^2)./(2*D) - 0.5;
% end