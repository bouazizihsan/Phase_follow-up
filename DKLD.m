function f=DKLD(u1,u2,sig1,sig2,k)
f  = 0.5 *(trace(inv(sig2)*sig1) + (u1 -u2)'*inv(sig2)*(u1-u2)  -k +log(det(sig2)/det(sig1)));
end