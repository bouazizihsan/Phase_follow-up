function u=maximum(Q)
    u=0;
    L=size(Q,1);
    v=Q(1,1);
    m=[];
    for i=1:L
        if(v<Q(i,2))
        v=Q(i,2);
        u=1;
        m=[m (i-1)];
        end
       
    end
    
end