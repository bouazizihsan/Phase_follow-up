% pour le QuatreQam
function   [xDec,xhat] = dem(y,M)
L = length(y);
xhat = zeros(1,L);
%bhat = [];
for i=1:L
    metric = zeros(1,M);
    for j = 1:M
        metric(j) =  norm( y(i) - mapping(j-1,M)).^2;
    end
    indxMin=find(metric==min(metric));
    xhat(i) = mapping(indxMin-1,M);
    xDec(i) = indxMin-1;
   % bhat = [bhat,de2bi(indxMin-1,2,'left-msb')];
end
end


function x = mapping(be,M)
   x = qammod(be,M);
end
function sig=signe(T)

sig=1;
for cc=1:length(T)
    sig=sig*sign(T(cc));
end

end
function som=somme(V)
som=0;
for k=1:length(V)
    som=som+V(k);
end
end
