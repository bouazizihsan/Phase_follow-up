function [Mnl,Sigma,PHI]=ProductProjection(Mnl_p0,Snl,Xci_p0,parameters1,parameters2,parameters3,InfP,AlphaGammaProj)
NClass=2;
M=size(Mnl_p0,1);
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

switch AlphaGammaProj
    case 'None'
        Mnl=m;
        Sigma=v;
        PHI=log(ph);
    case 'Classtring'
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
            Sigma=v;
            PHI=log(ph);
            
        else
            PHI=log(ph);
            Sigma=v;
            Mnl=m;
            
        end
    case 'Projection'
         [Mnl,Sigma]=Somme_Gaussian(m,v,ph);
%        [Mnl,Sigma]=Gaussian1(m,v,ph);
        PHI=0;
    otherwise  % case 'OneSideCorr' or 'BothSideCorr'
        if(M>1)
            ph=remove_irrelevant_terms( m, ph, InfP(1), InfP(2));
            [Mnl,Sigma]=Somme_Gaussian(m,v,ph./sum(ph));
            
            %             %         [Mnl,Sigma]=Gaussian1(m,v,ph);
            %             %             [Mnl,Sigma]=Somme_Gaussian(m,v,ph);
            %
%                         % ******** if we use F*alpha and Beta(Pilot) to find the best m= argmin(DKL(m,Mnl,v,Sigma) + DKL(m,InfP(1),v,InfP(2))) *********
%                         [Mnl,Sigma]=Somme_Gaussian(m,v,ph);
%                         f=(DKL(m,Mnl,v,Sigma) + DKL(m,InfP(1),v,InfP(2)));
%                         %             [val,ind]=min(f);
%                         %             Mnl=m(ind);
%                         %             Sigma=v(ind);
%                         [val ind] = sort(f);
%                         [Mnl,Sigma]=Somme_Gaussian(m(ind(1:NClass)),v(ind(1:NClass)),ph(ind(1:NClass))./sum(ph(ind(1:NClass))) );
            
            %             %                 [val,ind]=max(ph);
            %             %                 Sigma=v(ind);
            %             %                 Mnl=m(ind);
            %
            %             %                 % ******** if we use F*alpha*Beta(Pilot) *********
            %             %                 [val ind] = sort(ph);
            %             %                 [Mnl,Sigma]=Somme_Gaussian(m(ind(1:2)),v(ind(1:2)),ph(ind(1:2))./sum(ph(ind(1:2))));
            %
            %             %                 % ********END if we use F*alpha*Beta(Pilot) *********
            %
            %             %                 % ******** if we use F*alpha + Beta(Pilot) *********
            %             %                 [Mnl1,Sigma1]=Somme_Gaussian(m,v,ph);
            %             %                 [Mnl,Sigma]=Somme_Gaussian([Mnl1;InfP(1)],[Sigma1;InfP(2)],[0.5;0.5]);
            %             %                 Sigma=Sigma.* sqrt(InfP(2).*Sigma1);
            %             %                 % ********END if we use F*alpha+Beta(Pilot) *********
            PHI=0;
        else
            
            Mnl=m;
            Sigma=v;
            PHI=0;
        end
        
end

end
function f=DKL(A,B,C,D)
% dist=abs(Operator(B-A));
dist=abs((B-A));
f=dist.^2./(2*C) + D./(2*C) +0.5*log(C./D) - 0.5;
end

function phi_tilde=remove_irrelevant_terms( m_tilde, phi_tilde, mean_pilot_nn, var_pilot_nn )

pdf_Gauss   = @(y,m,s2) exp( -Operator(y-m).^2./(2*s2) )./sqrt(2*pi*s2); % pdf(y|p,x) real Gaussian
%pdf_Gauss   = @(y,m,s2) exp( -(y-m).^2./(2*s2) )./sqrt(2*pi*s2); % pdf(y|p,x) real Gaussian

CUT_threshold = 0.95; %% used if AlphaBetaPilotFilter='Y' to remove unnecessary elements from gamma
%PICK_f=0;

%% remove the irrelevant terms
phi_tilde_tmp = phi_tilde .* pdf_Gauss( m_tilde, mean_pilot_nn, var_pilot_nn  );

[vv, ii] = sort( phi_tilde_tmp/sum(phi_tilde_tmp) , 'descend');
vv_cum = cumsum(vv);
[fi]=find( vv_cum > CUT_threshold );
ii_zero = ii( fi(2:end) );
phi_tilde( ii_zero ) = 0;

%             switch PICK_f
%                 case 0%% the terms contributing to the distribution survive
%                     [fi]=find( vv_cum > CUT_threshold );
%                     ii_zero = ii( fi(2:end) );
%                 case 1%% only the maximum term survives
%                     ii_zero = ii( 2:end );
%                 case 2%%  the terms with the smallest mean survive
%                     [vv, ii]=sort( abs(m_tilde-mean_pilot_nn), 'ascend' );
%                     ii_zero = ii( 3:end );
%             end
end