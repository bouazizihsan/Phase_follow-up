function [MP,llrMP]= MsgPassing(LLR_codage,H,Imax)
metric=zeros(size(H,2),size(H,1));
r = metric;
A=100;
for comp=1:size(LLR_codage,1)% LLR est une matrice de taille: exp: 50*3 ; 50 fois la boucle ; 3 llr de code composé de 3 bits
    
    
        %initialisation
        for i=1:size(H,2)
            for j=1:size(H,1)
                if(H(j,i))
                    metric(i,j)=LLR_codage(comp,i);
                else
                    metric(i,j)=A;
                end
            end
        end
        
        %Horizontal step:check node j --->variable node i: check node apdate; r_ij ; dc: noeuds de check; dv: noeuds
        %des variables
        for n=1:Imax
%         r=zeros(size(H,2),size(H,1));
        for dc=1:size(H,1)% parcourir de check node
            
            for dv=1:size(H,2)% parcourir des variables node
                T = [];
                 if(H(dc,dv))
                for dvv=1:size(H,2)%parcourir les variables node sauf dv
                    if(dvv == dv)
                        T = [T,A];
                        continue;
                    end
                   if(H(dc,dvv))
                      T= [T,metric(dvv,dc)];
                   end

                end
                
                a=min(abs(T));
                sig=signe(T);
                r(dv,dc)=sig*a;
                 end
            end
        end
        % Vertical step: var node dv----> check node dc
        metric=zeros(size(H,2),size(H,1));
        for dv=1:size(H,2)
            for dc=1:size(H,1)
                V=[];
                V=[V,LLR_codage(comp,dv)];
                if(H(dc,dv))
                for dcc=1:size(H,1)%parcourir les
                    
                    if(dcc == dc)
                         V=[V,0];
                        
                        continue;     
                    end
                  
                        V=[V,r(dv,dcc)];
                    end
                
                metric(dv,dc)=somme(V);
                else
                    metric(dv,dc)=A;
                end
                
            end
        end
%         ___________________

        for dv=1:size(H,2)
            V=[];
            for dc=1:size(H,1)
                 V=[V,r(dv,dc)];
            end
                 V=[V,LLR_codage(comp,dv)];
                llrMP(comp,dv)=somme(V);
            xbit(dv)=1*(llrMP(comp,dv)>0)+0*(llrMP(comp,dv)<0);    
        end
         if(xbit*H'==0)
            break;   
        end   
        end
        
          MP(comp,:)=xbit;
    end

end



% _________________________________________________________________________

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