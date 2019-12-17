function jag=jac_log_Nclass2(x,coef,L,M)
for i=1:(M/2)
        z(:,i)=sum(coef.*exp(x(:,(i-1)*L+1:i*L )) ,2);
   
end
s=size(z);
y=max(z.').';
yy=repmat(y,1,s(2));
jag= y + log(sum(exp(-abs(z-yy)),2));
end