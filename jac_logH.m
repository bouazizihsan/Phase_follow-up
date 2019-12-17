function jag=jac_logH(x)
s=size(x);
y=max(x.').';
yy=repmat(y,1,s(2));
jag= y + log(sum(exp(-abs(x-yy)),2));
end