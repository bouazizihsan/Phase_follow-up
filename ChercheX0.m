function [x0,ind]= ChercheX0(g1,g2,m1,m2,phi1,phi2)
mx=max(m1,m2);
mi=min(m1,m2);

a=(g2-g1);

b=2*(m2.*g1 - m1.*g2);
c=(g2.*m1.^2-g1.*m2.^2 - 2*g1.*g2.*log(phi1./phi2));

delta=b.^2 - 4*a.*c;
ind=find(delta>0);
x0= (-b(ind)+sqrt(delta(ind) ))./(2*a(ind));

indx=find([x0<mi(ind)] + [x0>mx(ind)] > 0);
x0(indx)= (-b(indx)-sqrt(delta(indx)))./(2*a(indx));
i=find(a(ind)==0);
x0(i)=(2*g2(i).^2.*log(phi1./phi2(i))-(m1.^2 - m2(i).^2))./(2.*(m2(i)-m1));
end