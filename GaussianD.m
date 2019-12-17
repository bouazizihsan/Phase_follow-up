function [U_min, sigma2_min]=GaussianD(m,sigma1,phi)
s=size(m);
D=s(2);

Moyenne=[];
Variance=[];
DET=[];



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

for kk=1:s(1)+1
    
    
%     if(kk==(s(1)+2))
%         [mmm,fs]= sort(m(:,1));
%         phi2=phi(fs);
%         phi2=[phi2;phi2];
%         
%         ind=[D*(fs-1)+1 D*fs]';
%         ind=ind(:);
%         sigg=sigma1(ind,:);
%         sigg=[sigg;sigg];
%         Matx=[m(fs,:);m(fs,:)+2*pi];
%         %         Matx(:,2)=[Matx(1:s(1),2);Matx(s(1)+1:end,2)+2*pi];
%     else
        if(kk>1)
           fs2=[fs2(2:end);fs2(1)];        
        end
         [mmm,fsi]=sort(fs2);
          M2=vect2(kk:s(1)+kk-1,:);
          M2=M2(fsi);
          M2=M2(fs1);     
          Matx(:,2)=[M2;M2];
%     end
    
    
    
    for i=1:s(1)+1%i=(1+(kk>1)):s(1)
        
        moy=sum(repmat(phi2(i:s(1)+i-1),1,D).*Matx(i:s(1)+i-1,:),1)./(sum(repmat(phi2(i:s(1)+i-1),1,D),1));
        
        
        %         v1=(1:D:D*s(1)).';
        %         v2=(D:D:D*s(1)).';
        variable=zeros(D,D);
        mm=Matx(i:s(1)+i-1,:);
        
        q=phi2(i:s(1)+i-1);
        sig=sigg(1+(i-1)*D:s(1)*D+(i-1)*D,:);
        for j=1:s(1)
            variable=variable + q(j).*(mm(j,:).'*mm(j,:) + sig(1+(j-1)*D:j*D,:))./(sum(q));
        end
        var=variable - (moy.')*moy;
%         var=diag(diag(var));
        
        %         var(1+(i-1)*D:i*D,:)=([sum(variable(v1,:),1);sum(variable(v2,:),1)])./(sum(q)) - moy(i,:)*(moy(i,:)');
        det_min=(det(var));
        Moyenne=[Moyenne;moy];
        Variance=[Variance;var];
        DET=[DET;det_min];
    end
    
    
    
end


[mini,indx]=min(DET);
U_min=Moyenne(indx,:);
sigma2_min=Variance(1+(indx-1)*D:indx*D,:);

U_min(U_min>pi)=U_min(U_min>pi)-2*pi;
U_min(U_min<-pi)=U_min(U_min<-pi)+2*pi;

end
