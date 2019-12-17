function [cor,Mnl,S,PHI]=Projection_of_productD(Mnl_p0,Snl,Xci_p0,parameters1,parameters2,parameters3,Corr,AlphaGammaProj,D)

pdf_Gauss=@(t,m,s2,d) exp(-0.5*sum(((t-m )'*inv(s2).*(t-m ).'),2) )./((2*pi).^(d/2)*det(s2).^(0.5));
NClass=2;
epsilon=0.01;

PHI_nl_len_p0=[];
Mnl_tilde_p0=[];
Snl_tilde=[];
MM=size(Mnl_p0,1);
for j=1:size(parameters3,1)
    if(diag(parameters2((j-1)*D+1:j*D,:))==Inf)%(diag(alpha_p0{n-1,2})==(ones(D,1)*Inf))
        
        PHI_nl_len_p0=[PHI_nl_len_p0;(Xci_p0 + parameters3(j))];
        Mnl_tilde_p0 =[Mnl_tilde_p0;Mnl_p0];
        Snl_tilde    = [Snl_tilde; Snl];
        if(diag(Snl((MM-1)*D+1:MM*D,:))==Inf)
            PHI_nl_len_p0(end)=[];
            Mnl_tilde_p0(end,:)=[];
            Snl_tilde(end-D+1:end,:)= [];
        end
    else
        
        
        for i=1:MM
            [m,v,ph]=Produit_Alpha_Gamma_D(Mnl_p0(i,:) , Snl((i-1)*D+1:i*D,:) ,Xci_p0(i,:),parameters1(j,:),parameters2((j-1)*D+1:j*D,:),parameters3(j),D);
            PHI_nl_len_p0=[PHI_nl_len_p0;ph];
            Mnl_tilde_p0=[Mnl_tilde_p0;m];
            Snl_tilde=[Snl_tilde;v];
        end
        
        
        %% Case division
        %                 [m,v,ph]=Division_Alpha_Gamma_D(Mnl_p0 , Snl ,Xci_p0,parameters1,parameters2,parameters3,D);
        %                 PHI_nl_len_p0=[PHI_nl_len_p0;ph];
        %                 Mnl_tilde_p0=[Mnl_tilde_p0;m];
        %                 Snl_tilde=[Snl_tilde;v];
        
        
    end
end
PHI_nl_len_p0=exp(PHI_nl_len_p0(:,1)-jac_logH(PHI_nl_len_p0(:,1).'));
ind=find(PHI_nl_len_p0==0)';
if(ind)
    q=[(ind-1)*D+1;ind*D];
    q=q(:);
    
    PHI_nl_len_p0(ind,:)=[];
    Mnl_tilde_p0(ind,:)=[];
    Snl_tilde(q,:)=[];
    
end

% Mnl_tilde_p0(Mnl_tilde_p0<-pi)=Mnl_tilde_p0(Mnl_tilde_p0<-pi)+2*pi;
% Mnl_tilde_p0(Mnl_tilde_p0>pi)=Mnl_tilde_p0(Mnl_tilde_p0>pi)-2*pi;

switch AlphaGammaProj
    case 'None'
        PHI=log(PHI_nl_len_p0);
        S=Snl_tilde;
        Mnl=Mnl_tilde_p0;
        cor=1;
    case 'Classtring'
        %         PHI_nl_len_p0=exp(PHI_nl_len_p0(:,1)-jac_logH(PHI_nl_len_p0(:,1).'));
        % %         %         [val,ind]=max(PHI_nl_len_p0);
        % %         %         q=[(ind-1)*D+1 ind*D]';
        % %         %         q=q(:);
        % %         %         S=Snl_tilde(q,:);
        % %         %         Mnl=Mnl_tilde_p0(ind,:);
        % %         %         PHI=1;
        
        %_______________________Max__________________________
        %         if(length(PHI_nl_len_p0)>1)
        %             [val ind] = sort(PHI_nl_len_p0,'descend');
        %             q=[(ind(1:NClass)-1)*D+1 ind(1:NClass)*D]';
        %             q=q(:);
        %             S=Snl_tilde(q,:);
        %             Mnl=Mnl_tilde_p0(ind(1:NClass),:);
        %             PHI=PHI_nl_len_p0(ind(1:NClass));
        %             PHI=PHI./sum(PHI);
        %             PHI=log(PHI);
        %         else
        %             S=Snl_tilde;
        %             Mnl=Mnl_tilde_p0;
        %             PHI=log(PHI_nl_len_p0);
        %         end
        %         cor=1;
        MM=length(PHI_nl_len_p0);
        if(MM>1)
            f=[];
            indx1=[];
            
            while(length(PHI_nl_len_p0)>NClass)
                if(norm(f))
                    indice= [repmat(length(PHI_nl_len_p0),length(PHI_nl_len_p0)-1,1),(1:(length(PHI_nl_len_p0)-1))'];
                else
                    indice=(1:length(PHI_nl_len_p0));
                    indice=nchoosek(indice,2);
                end
                Mnl=zeros(length(indice),D);
                S=zeros(length(indice)*D,D);
                for j=1:length(indice)
                    
                    q=[(indice(j,:)-1)*D+1;indice(j,:)*D];
                    q=q(:);
                    [Mnl(j,:),S((j-1)*D+1:j*D,:)]=Somme_GaussianD(Mnl_tilde_p0(indice(j,:),:),Snl_tilde(q,:),PHI_nl_len_p0(indice(j,:)));
                    
                    f=[f;sum(PHI_nl_len_p0(indice(j,:))./(sum(PHI_nl_len_p0(indice(j,:)))).*DKLD(Mnl_tilde_p0(indice(j,:),:).',Mnl(j,:).',Snl_tilde(q,:),S((j-1)*D+1:j*D,:),D))];
                    if( isnan(norm(f)) )
                        stp=1;
                    end
                end
                
                if(length(indx1))
                    indice=[indx1;indice];
                    S=[S1;S];
                    Mnl=[Mnl1;Mnl];
                end
                [val ind]=min(f);
                q=(indice(ind,:)).';
                [X,Y]=find((indice==q(1)) + (indice==q(2)));
                f(X)=[];
                
                coef=sum(PHI_nl_len_p0(indice(ind,:)));
                q=[(indice(ind,:)-1)*D+1;indice(ind,:)*D];
                q=q(:);
                Snl_tilde(q,:)=[];
                Mnl_tilde_p0(indice(ind,:),:)=[];
                PHI_nl_len_p0(indice(ind,:))=[];
                Snl_tilde=[Snl_tilde;S((ind-1)*D+1:ind*D,:)];
                Mnl_tilde_p0=[Mnl_tilde_p0;Mnl(ind,:)];
                PHI_nl_len_p0=[PHI_nl_len_p0;coef];
                S((X-1)*D+1:X*D,:)=[];
                S1=S;
                Mnl(X,:)=[];
                Mnl1=Mnl;
                if(length(PHI_nl_len_p0)>2)
                    indx1=(1:(length(PHI_nl_len_p0)-1)).';
                    indx1=nchoosek(indx1,2);
                else
                    indx1=[];
                end
            end
            Mnl=Mnl_tilde_p0;
            S=Snl_tilde;
            PHI=log(PHI_nl_len_p0);
            cor=1;
        else
            PHI=log(PHI_nl_len_p0);
            S=Snl_tilde;
            Mnl=Mnl_tilde_p0;
            cor=1;
        end
        
        
        %______________________________KLTRack______________________
        %                 j=1;
        %                 S=[];
        %                 Mnl=[];
        %                 PHI=[];
        %                 while ((j <= NClass) && (norm(PHI_nl_len_p0)>0))
        %                     [var, indx] = max(PHI_nl_len_p0);
        %
        %                     f=DKLD(Mnl_tilde_p0.',Mnl_tilde_p0(indx,:).',Snl_tilde,Snl_tilde((indx-1)*D+1:indx*D,:),D);
        %                     aux=find(f<epsilon);
        %                     sum_phi=sum(PHI_nl_len_p0(aux));
        %
        %         %             q=[(indx-1)*D+1 indx*D]';
        %         %             q=q(:);
        %         %             S= [S;Snl_tilde(q,:)];
        %         %             Mnl= [Mnl;Mnl_tilde_p0(indx,:)];
        %         %             PHI= [PHI;sum_phi];
        %
        %                     q=[(aux-1)*D+1 aux*D]';
        %                     q=q(:);
        %                     [auxM,auxS]=GaussianD1(Mnl_tilde_p0(aux,:),Snl_tilde(q,:),PHI_nl_len_p0(aux));
        %                     S= [S;auxS];
        %                     Mnl= [Mnl;auxM];
        %                     PHI= [PHI;sum_phi];
        %
        %                     Mnl_tilde_p0(aux,:)= [];
        %                     Snl_tilde((aux-1)*D+1:aux*D,:)=[];
        %                     PHI_nl_len_p0(aux)=[];
        %                     j=j+1;
        %                 end
        %                  cor=sum(PHI);
        %                 if((1-sum(PHI))>(10.^(-20))) % 10.^(-323) limit de calcul de log
        %         %             [auxM,auxS]=GaussianD1(Mnl_tilde_p0,Snl_tilde,PHI_nl_len_p0);
        %         %             S= [S;auxS];
        %         %             Mnl= [Mnl;auxM];
        %         %             PHI= [PHI;sum(PHI_nl_len_p0)];
        %
        %                     S= [S;diag(Inf*(1:D))];
        %                     Mnl= [Mnl;0.*(1:D)];
        %                     PHI= [PHI;(1-sum(PHI))];
        %                 end
        %
        
        % %         if(norm(PHI_nl_len_p0))
        % %
        % %             [auxM,auxS]=GaussianD1(Mnl_tilde_p0,Snl_tilde,PHI_nl_len_p0);
        % %
        % %             sum_phi= sum(PHI_nl_len_p0);
        % %             %         if(S<3)
        % %             S= [S;auxS];
        % %             Mnl= [Mnl;auxM];
        % %             PHI= [PHI;sum_phi];
        % %             %         else
        % %             %             S= [S;Inf];
        % %             %         Mnl= [Mnl,auxM];
        % %             %         PHI= [PHI;sum_phi];
        % %             %         end
        % %             indx=1;
        % %         end
        %
        %         PHI=PHI./sum(PHI);
        %
        %         PHI=log(PHI);
        
        % %         HB-Track
        %                 [ind val] = find(PHI_nl_len_p0>epsilon);
        %
        %                 Mnl_tilde_p0=Mnl_tilde_p0(ind,:);
        %                 q=[(ind-1)*D+1 ind*D]';
        %                 q=q(:);
        %                 Snl_tilde=Snl_tilde(q,:);
        %                 PHI_nl_len_p0=PHI_nl_len_p0(ind,:);
        %                 cor=sum(PHI_nl_len_p0);
        %                 PHI=PHI_nl_len_p0./sum(PHI_nl_len_p0);
        %                 [Mnl,S]=GaussianD1(Mnl_tilde_p0,Snl_tilde,PHI_nl_len_p0);
        %                 PHI=0; %log(1)
        
    case 'Projection'
        [Mnl,S]= Somme_GaussianD(Mnl_tilde_p0,Snl_tilde,PHI_nl_len_p0) ;
        cor=1;
        PHI=0;
    otherwise  % case 'OneSideCorr' or 'BothSideCorr'
        MM=size(Mnl_tilde_p0,1);
        if(MM>1)
            %% By DKL
%             f=[];
%             %                 PHI_nl_len_p0=exp(PHI_nl_len_p0(:,1)-jac_logH(PHI_nl_len_p0(:,1).'));
%             [Mnl,S]= Somme_GaussianD(Mnl_tilde_p0,Snl_tilde,PHI_nl_len_p0) ;
%             
%             for j=1:MM
%                 q=[(j-1)*D+1;j*D];
%                 q=q(:);
%                 f=[f;(DKLD(Mnl_tilde_p0(j,:).',Mnl.',Snl_tilde(q,:),S,D)+DKLD(Mnl_tilde_p0(j,:).',Corr{1,1}.',Snl_tilde(q,:),Corr{2,1},D))];
%             end
%             [val ind] = sort(f);
%             q=[(ind(1:NClass)-1)*D+1 ind(1:NClass)*D]';
%             q=q(:);
%             [Mnl,S]= Somme_GaussianD(Mnl_tilde_p0(ind(1:NClass),:),Snl_tilde(q,:),PHI_nl_len_p0(ind(1:NClass),:)./repmat(sum(PHI_nl_len_p0(ind(1:NClass),:)),NClass,1)) ;
            %% By remove_irrelevant_terms
            PHI_nl_len_p0=remove_irrelevant_termsD( Mnl_tilde_p0, PHI_nl_len_p0, Corr{1,1}, Corr{2,1});
            [Mnl,S]=Somme_GaussianD(Mnl_tilde_p0,Snl_tilde,PHI_nl_len_p0) ;
        else
            Mnl=Mnl_tilde_p0;
            S=Snl_tilde;
        end
        cor=1;
        PHI=0;
end

end

% function f=DKL(A,B,C,D)
% %f= sum(A.*log(A./repmat(B,1,size(A,2))),1);
% dist=abs(Operator(B-A));
% 
% f=dist.^2./C + D./C +0.5*log(C./D) - 1;
% 
% end
function f=DKLD(u1,u2,sig1,sig2,D)
f=[];
for i=1:size(u1,2)
    f=[f; (0.5 *(trace(inv(sig2)*sig1((i-1)*D+1:i*D,:)) + (u1(:,i) - u2)'*inv(sig2)*(u1(:,i) - u2)  - D + log(det(sig2)/det(sig1((i-1)*D+1:i*D,:)))))];
end
end

function phi_tilde=remove_irrelevant_termsD( m_tilde, phi_tilde, mean_pilot_nn, var_pilot_nn )

pdf_Gauss=@(t,m,s2,d) exp(sum((-0.5*((t-m ))'*inv(s2).*((t-m ))),2) )./( (2*pi).^(d/2)*det(s2).^(0.5) );

CUT_threshold = 0.95; %% used if AlphaBetaPilotFilter='Y' to remove unnecessary elements from gamma
%PICK_f=0;
D=2;
M=size(m_tilde,1);
%% remove the irrelevant terms
phi_tilde_tmp = phi_tilde .* pdf_Gauss( m_tilde.', repmat(mean_pilot_nn.',1,M), var_pilot_nn,D);

[vv, ii] = sort( phi_tilde_tmp/sum(phi_tilde_tmp) , 'descend');
vv_cum = cumsum(vv);
[fi]=find( vv_cum > CUT_threshold );
ii_zero = ii( fi(2:end) );
phi_tilde( ii_zero ) = 0;
phi_tilde=phi_tilde./sum(phi_tilde);
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
