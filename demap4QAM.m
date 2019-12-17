% Hsan; =dem
function vectb = demap4QAM(vects)
l = length(vects); 
jj=1;
ll=l.*4;
for i=1:l
    
    mut(jj)=(real(vects(i))-1).^2 + (imag(vects(i))-1).^2;
    mut(jj+1)=(real(vects(i))+1).^2 + (imag(vects(i))-1).^2;
    mut(jj+2)=(real(vects(i))-1).^2 + (imag(vects(i))+1).^2;
    mut(jj+3)=(real(vects(i))+1).^2 + (imag(vects(i))+1).^2;
    jj=jj+4;
end
minarg = min(mut);
kk=1;
for ii=1:length(minarg); 
    switch minarg(ii)
        case 0
            vectb(kk)=0;
            vectb(kk+1)=0;
            kk=kk+2;
        case 1
            vectb(kk)=0;
            vectb(kk+1)=1;
            kk=kk+2;
        case 2
            vectb(kk)=1;
            vectb(kk+1)=0;
            kk=kk+2;
        case 3
            vectb(kk)=1;
            vectb(kk+1)=1;
            kk=kk+2;
    end
end

end 


function  minarg =min(mut)
k=1;
ll = length(mut);
for i=1:4:ll
    
minarg(k)=0;

for j=i+1:4
    if(mut(j)<mut(j-1))
        minarg(k)=j-1; 
    end
 k=k+1; 
end 
end
end