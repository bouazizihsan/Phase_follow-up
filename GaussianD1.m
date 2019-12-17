function [U_min, sigma2_min]=GaussianD1(m,sigma1,phi)
s=size(m);
D=s(2);
phi=phi./sum(phi);
Moyenne=[];
Variance=[];
DET=[];
dig=[];


[mmm,fs1]= sort(m(:,1));
phi2=phi(fs1);
phi2=[phi2;phi2];

ind=[D*(fs1-1)+1 D*fs1]';
ind=ind(:);
sigg=sigma1(ind,:);
sigg=[sigg;sigg];
Matx=[m(fs1,:);m(fs1,:)];
Matx(:,1)=[m(fs1,1);m(fs1,1)+2*pi];

[mmm,fs2]= sort(m(:,2));
vect2=[m(fs2,2);m(fs2,2)+2*pi];

for kk=1:s(1)
        if(kk>1)
           fs2=[fs2(2:end);fs2(1)];        
        end
         [mmm,fsi]=sort(fs2);
          M2=vect2(kk:s(1)+kk-1,:);
          M2=M2(fsi);
          M2=M2(fs1);     
          Matx(:,2)=[M2;M2];
    
    
    for i=1:s(1)
        
        moy=sum(repmat(phi2(i:s(1)+i-1),1,D).*Matx(i:s(1)+i-1,:),1)./(sum(repmat(phi2(i:s(1)+i-1),1,D),1));

        variable=zeros(D,D);
        mm=Matx(i:s(1)+i-1,:);
        
        q=phi2(i:s(1)+i-1);
        sig=sigg(1+(i-1)*D:s(1)*D+(i-1)*D,:);
        
        for j=1:s(1)
            variable=variable + q(j).*((mm(j,:)-moy).'*(mm(j,:)-moy) + sig(1+(j-1)*D:j*D,:));%./(sum(q));
        end
        var=variable;% - (moy.')*moy;
%         var=diag(diag(var));
        det_min=(det(var));
        dig=[dig;prod(diag(var))];
        Moyenne=[Moyenne;moy];
        Variance=[Variance;var];
        DET=[DET;det_min];
    end
    
    
    
end

% [mini,indx]=min(DET);
[mini,indx]=min(dig);
if((indx~=1) && (indx~=4))
   stp=1; 
end
U_min=Moyenne(indx,:);
sigma2_min=Variance(1+(indx-1)*D:indx*D,:);

U_min(U_min>pi)=U_min(U_min>pi)-2*pi;
U_min(U_min<-pi)=U_min(U_min<-pi)+2*pi;

end
