function jag=jac_log_Nclass(x,coef,coef_INF,L,M)
for i=1:(M/2)
    if(norm(coef_INF))
        z(:,i)=sum(coef.*exp(x(:,(i-1)*L+1:i*L )) + sum(coef_INF,2)./(4*pi^2) ,2);
    else
        z(:,i)=sum(coef.*exp(x(:,(i-1)*L+1:i*L )) ,2);
    end
end
s=size(z);
y=max(z.').';
yy=repmat(y,1,s(2));
jag= y + log(sum(exp(-abs(z-yy)),2));
end