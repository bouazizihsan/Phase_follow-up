function [betaj,betajG,zT,zG]=KLCluster_U_new(f,f2,z_T,z_G,epsilon,Len)

phi=f(:);
phi=phi./sum(phi);
z=z_T(:);
z_GG=z_G(:);
j=1;
betaj=[];
betajG=[];
zT=[];
zG=[];


while((j<=Len) && (norm(phi)>0) )
    [lead,ind]=max(phi);
    dist=DKL(z(ind),z);
    aux=find(dist<epsilon);
    betaj=[betaj;sum(phi(aux))];
    zz=z(aux);
    var=phi(aux)./betaj(j);
    [zn]=CMVM(var,zz,1);
    zT=[zT;zn];
    phi(aux)=[];
    z(aux)=[];
    j=j+1;
end
fact= 1 - (sum(betaj));

if((norm(phi)) && (fact>0))
        phi=phi./sum(phi);
        betaj=[betaj;fact];
        [zn]=CMVM(phi,z,1);
        zT=[zT;zn];
 
end

%_____________

phi=f2(:);
phi=phi./sum(phi);
j=1;
while((j<=Len) && (norm(phi)>0) )
    [lead,ind]=max(phi);
    dist=DKL(z_GG(ind),z_GG);    
    aux=find(dist<epsilon);
    betajG=[betajG;sum(phi(aux))];
    zz=z_GG(aux);
    var=phi(aux)./betajG(j);
    [zn]=CMVM(var,zz,2);
   
    zG=[zG;zn];
    phi(aux)=[];
    z_GG(aux)=[];
    j=j+1;
end
fact= 1 - (sum(betajG));

if((norm(phi)) && (fact>0))
        phi=phi./sum(phi);
        betajG=[betajG;fact];
        [zn]=CMVM(phi,z_GG,2);
        
        zG=[zG;zn];   
end

%___________________

end

function [g]=CMVM(alph,z,i)

alpha=alph/sum(alph);
if(i==1)
u=angle(sum(alpha.*(1-1./(2*abs(z) ) ).*exp(1j*angle(z)) ) );
k=1./(sum(alpha./abs(z)));
g=k.*exp(1j.*u);
else
m=angle(z);
var1=1-(besseli(1,abs(z))./besseli(0,abs(z)));
%moy=sum(m.*alpha/sum(alpha));
%var=sum( (alpha/sum(alpha)) .* ( var1 + m.^2)) - moy.^2;
[moy, var]=Gaussian1(m,var1,alpha);
kk= 1./(2*var);
uu=angle(besseli(1,kk)./ besseli(0,kk) .* exp(1j.*moy)) ; % =moy
% [uu moy u];
% [var 1./[kk k]];
g=kk .* exp(1j.*uu);
end

end

function d=DKL(A,B)
AA=abs(A);
BB=abs(B);
ua=angle(A);
ub=angle(B);
d=log(besseli(0,AA)./besseli(0,BB)) + besseli(1,BB)./besseli(0,BB).*(BB - AA.*cos(ua-ub));

end