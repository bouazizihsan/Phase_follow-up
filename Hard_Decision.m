function [X_Hard,i]=Hard_Decision(y,x)

M=length(x);
xx=repmat(x,length(y),1);
yy=repmat(y,1,M);

D=abs(yy-xx).^2; %% distance mesure
[v,i]=min(D.');
X_Hard=x(i);

end