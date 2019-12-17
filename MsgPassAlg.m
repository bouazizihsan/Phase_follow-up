function val_decision= MsgPassAlg(LLR,H,Imax)
% H=[1 1 1 0 1 0 0; 0 1 1 1 0 1 0; 1 1 0 1 0 0 1];
s=sparse(H);
[lines, columns, values] = find(s);
T=zeros(2,length(lines));
T(1,:)=lines';
T(2,:)=columns';
T=T(:)';
var_check=zeros(1,length(lines)); %tableau des metrics de variable vers check nodes
check_var=zeros(1,length(lines));%tableau des metrics check nodes vers les variables
%initialisation
for i=1:length(var_check)
    var_check(i)=LLR(columns(i));
end

%Horizontal step:check node j --->variable node i: check node apdate; r_ij ;

for n=1:Imax
    for j=1:length(check_var)
        % check_var(j)= check T(2*j-1)---->var T(2*j)
        %trouver l'occurence de check node dans  T ---> trouver les variables
        %liée à elle
        indx_check=[];
        indx_check=[indx_check,find(T==T(2*j-1))];
        indx=[];
        for c=1:length(indx_check)
            if(mod(indx_check(c),2)>0 && indx_check(c)~= (2*j-1) )
                indx= [indx,indx_check(c)+1];% indices des variables en relation avec check T(2*j-1)
            end
        end
        aux=[];
        for k=1:length(indx)
            aux=[aux,var_check(indx(k)./2)];
        end
        a=min(abs(aux));
        sig=signe(aux);
        
        check_var(j)=sig*a;
    end
    
    
    % Vertical step: var node ----> check node
    val_decision=[];
    for i=1:length(var_check)
        % var_check(j)= var T(2*j)---->check T(2*j-1)
        %trouver l'occurence de var node dans  T ---> trouver les check node
        %liée à elle
        indx_var=[];
        indx_var=[indx_var,find(T==T(2*i))];
        indx=[];
        indxtot=[];
        for c=1:length(indx_var)
            if((mod(indx_var(c),2)==0))
                
                if((indx_var(c)~= (2*i)))
                    indx= [indx,indx_var(c)];% indx_check(c)-1 :indices des check nodes en relation avec var T(2*j)
                end
                indxtot=[indxtot,indx_var(c)];
            end
            
        end
        aux=[];
        aux2=[];
        aux=[aux,LLR(columns(i))];
        aux2=[aux2,LLR(columns(i))];
        for k=1:length(indx)
            aux=[aux,check_var(indx(k)./2)];
        end
        for comp=1:length(indxtot)
            aux2=[aux2,check_var(indxtot(comp)./2)];
        end
        
        var_check(i)=somme(aux);
        if(i>1 && (columns(i)==columns(i-1)))
            continue ;
        end
        val_decision=[val_decision,somme(aux2)];
    end
    
end

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
