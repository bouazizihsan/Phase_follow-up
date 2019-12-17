function LLR=LLR_PN_D(Y, N0, BCJR, hMap1,hMap2,payload,TYPE, SUM )


if nargin<8, SUM='sumexp'; end

XX=[hMap1.Xsort hMap2.Xsort];
M=size(XX,1);
S=size(Y);
D=S(2);
k10=hMap1.pntk0;
k11=hMap1.pntk1;
% % k20=hMap2.pntk0;
% % k21=hMap2.pntk1;

[K1,kM]=size(hMap1.pntk0);
[K2,kM]=size(hMap2.pntk0);

if (K1~=M/2 || K2~=M/2) , error('number of rows in kk0 must be equal to M/2'); end
LLR=zeros(length(Y),kM);
NClass=2;
Theta_zero = BCJR.P_Ref;

LLR_log=zeros(S(1),M);
switch TYPE
    case 'SP1'
        for j=1:S(1)
            variance_PN =BCJR.alpha_beta{j,2};
            %     moyenne_PN  =BCJR.alpha_beta{j,1};
            Poids=BCJR.alpha_beta{j,3};
            for i=1:M
                variable=0;
                for k=1:1  %NClass.^2
                    
                    %         P_Ref= angle(y./x)
                    Diff=( (Theta_zero(i,:,j)).'); %-moyenne_PN(k,:)).');
                    %             Diff(Diff>pi)=Diff(Diff>pi)-2*pi;
                    %             Diff(Diff<-pi)=Diff(Diff<-pi)+2*pi;
                    variable(k)=((Poids(k)) - (abs(Y(j,:))-abs(XX(i,:)))*inv(N0)*(abs(Y(j,:))-abs(XX(i,:))).'...
                        - (0.5)*log(det((N0+2*diag(abs(Y(j,:)).*abs(XX(i,:)))*variance_PN((k-1)*S(2)+1:k*S(2),:) ) )) ...
                        + (-0.5*((Diff)'*2*diag(abs(Y(j,:)).*abs(XX(i,:)))*inv(2*diag(abs(Y(j,:)).*abs(XX(i,:)))*variance_PN((k-1)*S(2)+1:k*S(2),:) + N0)*Diff)) );
                    
                end
                LLR_log(j,i)=variable;%jac_log(variable); %pour k=1:NClass^2
                
            end
        end
    case 'LT'
        N0=diag(N0)';
        for j=1:S(1)
            variance_PN =BCJR.alpha_beta{j,2};
            %     moyenne_PN  =BCJR.alpha_beta{j,1};
            Poids=BCJR.alpha_beta{j,3};
            for i=1:M
                variable=0;
                for k=1:1  %NClass.^2
                    %                     variable(k)=(Poids(k)) - (abs(Y(j,:))-abs(XX(i,:)))*inv(N0)*(abs(Y(j,:))-abs(XX(i,:))).'...
                    %                         - (0.5)*log(det( (N0 + 2*diag(abs(XX(i,:)).^2)*variance_PN((k-1)*S(2)+1:k*S(2),:) ) )) ...
                    %                         + (0.5*(imag(conj(XX(i,:)).*Y(j,:)).')'*(2*variance_PN((k-1)*S(2)+1:k*S(2),:)*inv( 2*diag(abs(XX(i,:)).^2)*N0*variance_PN((k-1)*S(2)+1:k*S(2),:) + N0.^2 )*(imag(conj(XX(i,:)).*Y(j,:)).') ));
                    
                    variable(k)=(Poids(k))+ sum(- 2*log(abs(XX(i,:))) - 0.5*log(N0./abs(XX(i,:).^2)) - (real(Y(j,:)./XX(i,:) - 1)-1).^2 ./(N0./abs(XX(i,:).^2)))...
                        -0.5*log(det( diag(N0./(2*abs(XX(i,:).^2)))+ variance_PN((k-1)*S(2)+1:k*S(2),:) ))...
                        -0.5*(imag(Y(j,:)./XX(i,:)).')'*inv(diag(N0./(2*abs(XX(i,:).^2)))+ variance_PN((k-1)*S(2)+1:k*S(2),:) )*(imag(Y(j,:)./XX(i,:)).');
                end
                LLR_log(j,i)=variable;%jac_log(variable); %pour k=1:NClass^2
                
            end
        end
    otherwise
        error('unknown integral simplification');
end
% XX=XX.';
% yy=repmat(Y,1,1,M);% (Ns,D,M)
% Xref=repmat(XX,Ns,1,1);% (Ns,D,M)
% % variance_PN =BCJR.alpha_beta{:,2};
% moyenne_PN  =BCJR.alpha_beta{:,1};
% moyenne_PN  =repmat(moyenne_PN,1,1,M);
% P_Ref=BCJR.Mnl;
% P_Reff=repmat(P_Ref,1,1,M);
% %case 'SP1'
% y2 = conj(yy).*exp(1j*P_Reff); xref2 = xref;
% LLR_log=-0.5*sum(log(pi*N0.*(abs(xref2).^2))-abs(real(y2.*xref2)-(abs(xref2).^2)).^2./(N0.*(abs(xref2).^2)),2);  %% Euclidean distance
%
% delta =P_Reff - imag(y2./conj(xref));
% Diff=moyenne_PN-delta;
% Diff(Diff>pi)=Diff(Diff>pi)-2*pi;
% Diff(Diff<-pi)=Diff(Diff<-pi)+2*pi;
% vVar_PN =BCJR.var;
% vVar_PN =repmat(vVar_PN,1,1,1,M);
%
%     v2 = ( (N0.*(abs(xref2).^2)./2)+vVar_PN );
%     LLR_log=LLR_log-0.5*log(2*pi*v2)-(Diff).^2./(2.*v2) ;



for k=1:kM
    switch SUM
        case 'sumexp'
            LLR(:,k)=jac_log( LLR_log(:,k11(:,k)) ) - jac_log( LLR_log(:,k10(:,k)) );
        case 'maxlog'
            LLR(:,k)= max( LLR_log(:,k11(:,k)), [],2 ) - max( LLR_log(:,k10(:,k)), [],2);
        otherwise
            error('unknown metric combining');
    end
end
end
%%% en utilisant sigma= BCJR.Snl;
%             Diff=( (Theta_zero(i,:,j)-moyenne_PN(k,:)).');
%             Diff(Diff>pi)=Diff(Diff>pi)-2*pi;
%             Diff(Diff<-pi)=Diff(Diff<-pi)+2*pi;
%             variable(k)=((Poids(k)) - log(prod(abs(XX(i,:)).^2)) ...
%                 + (-0.5*(real(Y(j,:)./XX(i,:))-1)*inv(sigma{1,payload(j)}((i-1)*S(2)+1:i*S(2),:))*(real(Y(j,:)./XX(i,:)).'-1))...
%                 - (0.5)*log(det(2*pi*(sigma{1,payload(j)}((i-1)*S(2)+1:i*S(2),:)+variance_PN((k-1)*S(2)+1:k*S(2),:)))) ...
%                 + (-0.5*((Diff)'*inv(variance_PN((k-1)*S(2)+1:k*S(2),:) + sigma{1,payload(j)}((i-1)*S(2)+1:i*S(2),:))*Diff)) );


