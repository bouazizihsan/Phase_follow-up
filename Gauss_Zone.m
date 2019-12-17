function [U_min, sigma2_min]=Gauss_Zone(m,sigma1,phi)
s=size(m);
D=s(2);

Moyenne=[];
Variance=[];
DET=[];




for kk=1:D+2
    
%     if(kk==D+2)
%        [mmm,fs]= sort(m(:,2));
%         phi2=phi(fs);
%         phi2=[phi2;phi2];
%         
%         ind=[D*(fs-1)+1 D*fs]';
%         ind=ind(:);
%         sigg=sigma1(ind,:);
%         sigg=[sigg;sigg];
%         Matx(:,:,kk*D-2)=[m(fs,:);m(fs,:)+2*pi]; %Matx(:,:,kk*D-2)=[[m(fs,1);m(fs,1)+2*pi] [m(fs,2)+2*pi;m(fs,2)]];%6
%     else
   if(kk==D+1)
       [mmm,fs]= sort(m(:,1));
        phi2=phi(fs);
        phi2=[phi2;phi2];
        
        ind=[D*(fs-1)+1 D*fs]';
        ind=ind(:);
        sigg=sigma1(ind,:);
        sigg=[sigg;sigg];
        Matx(:,:,kk*D-1)=[m(fs,:);m(fs,:)+2*pi]; %5
%         Matx(:,2)=[Matx(1:s(1),2);Matx(s(1)+1:end,2)+2*pi];
    else
        [mmm,fs]= sort(m(:,kk));
        phi2=phi(fs);
        phi2=[phi2;phi2];
        
        ind=[D*(fs-1)+1 D*fs]';
        ind=ind(:);
        sigg=sigma1(ind,:);
        sigg=[sigg;sigg];
        Matx(:,:,kk*D-1)=[m(fs,:);m(fs,:)];
        Matx(:,kk,kk*D-1)=[m(fs,kk);m(fs,kk)+2*pi];
        Matx(:,:,kk*D)=[m(fs,:)+2*pi;m(fs,:)+2*pi];
        Matx(:,kk,kk*D)=[m(fs,kk);m(fs,kk)+2*pi];
    end
end

for kk=1:size(Matx,3)
    for i=1:s(1)%i=(1+(kk>1)):s(1)
        
        moy=sum(repmat(phi2(i:s(1)+i-1),1,D).*Matx(i:s(1)+i-1,:,kk),1)./(sum(repmat(phi2(i:s(1)+i-1),1,D),1));
        
        
        %         v1=(1:D:D*s(1)).';
        %         v2=(D:D:D*s(1)).';
        variable=zeros(D,D);
        mm=Matx(i:s(1)+i-1,:,kk);
        
        q=phi2(i:s(1)+i-1);
        sig=sigg(1+(i-1)*D:s(1)*D+(i-1)*D,:);
        for j=1:s(1)
            variable=variable + q(j).*(mm(j,:).'*mm(j,:) + sig(1+(j-1)*D:j*D,:))./(sum(q));
        end
        var=variable - (moy.')*moy;
         var=diag(diag(var));
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
