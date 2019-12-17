function [betaj,zT]=KLCluster_U(f,z_T,epsilon,Len,Threshold)

phi=f(:);
phi=phi./sum(phi);
z=z_T(:);
j=1;
betaj=[];
zT=[];

while((j<=Len) && (norm(phi)>0) )
    [lead,ind]=max(phi);
    dist=DKL(z(ind),z);
    aux=find(dist<epsilon);
    betaj=[betaj;sum(phi(aux))];
    zz=z(aux);
    var=phi(aux)./betaj(j);
    [zn]=CMVM(var,zz,Threshold);
    zT=[zT;zn];
    phi(aux)=[];
    z(aux)=[];
    j=j+1;
end

fact= 1 - (sum(betaj));

if((norm(phi)) && (fact>0))
        phi=phi./sum(phi);
        betaj=[betaj;fact];
        [zn]=CMVM(phi,z,Threshold);
        zT=[zT;zn];
 
end

end



function d=DKL(A,B)
AA=abs(A);
BB=abs(B);
ua=angle(A);
ub=angle(B);
d=log(besseli(0,AA)./besseli(0,BB)) + besseli(1,BB)./besseli(0,BB).*(BB - AA.*cos(ua-ub));

end