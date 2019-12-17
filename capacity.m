clc;
clear all;
warning off;
Ixy=0;
X=randi([0,1],1,2);% x(ii)--> Pr(i)
Y=randi([0,1],1,2);% y(jj)--> Pr(j) meme vecteur 
Pr=(0:0.1:1);
err=0.0625;
err2=0,25;
for i=1:length(Pr)
    i
    Ixy(i)=-(Pr(i)*(1-err)+(1-Pr(i))*err)*log2(Pr(i)*(1-err)+(1-Pr(i))*err)-(Pr(i)*err+(1-Pr(i))*(1-err))*log2((Pr(i)*err+(1-Pr(i))*(1-err)))+(1-err)*log2(1-err)+err*log2(err);  
    Ixy2(i)=-(Pr(i)*(1-err2)+(1-Pr(i))*err2)*log2(Pr(i)*(1-err2)+(1-Pr(i))*err2)-(Pr(i)*err2+(1-Pr(i))*(1-err2))*log2((Pr(i)*err2+(1-Pr(i))*(1-err2)))+(1-err2)*log2(1-err2)+err2*log2(err2);   
end
  
plot(Pr,Ixy) 
hold on
plot(Pr,Ixy2) 
    

