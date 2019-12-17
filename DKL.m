function f=DKL(A,B,C,D)
%f= sum(A.*log(A./repmat(B,1,size(A,2))),1);
dist=abs(B-A);
[ind v]=find(dist>pi);
dist(ind)=dist(ind)-2*pi;

f=dist.^2./(2*C) + D./(2*C) + 0.5*log(C./D) - 0.5;

end