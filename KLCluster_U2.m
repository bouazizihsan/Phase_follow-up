function [Coefs,zT]=KLCluster_U2(f,z_T,NClass,Threshold)
z=z_T(:);
phi=f(:);
PH=phi./sum(phi);
Moy=angle(z);
S = 1./(abs(z));

[zT,Coefs] = DKL_X0(Moy,S,PH,NClass,Threshold,'Tikhonov');

end
