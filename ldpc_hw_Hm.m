function Hobj = ldpc_hw_Hm(setup)
%% *********************************************************************
%%% The matrices are MxN sub-matrices of ZxZ elements each
%%% three original matrices OM all have the same number of columns N=96
%%% the number of row M depends on the rate
%%% R=1/2--> Hb12(M=48), R=3/4--> Hb34(M=24), R=7/8--> Hb78(M=12)
%%% each element of OM is a matrix of Per=ZxZ elements
%%% Z=84 for Blength=8064, Z=42 for Blength=4032, Z=21 for Blength=2016
%%% Per matrix is a permutation of identity matrix
%%% For each rate the same Original matrix for Blength = 8064 can be used
%%% other Blengths, the rule is HbNew(n,m)=floor(HbOrig(n,m)/zoom); zoom=2
%%% for Z=42, and zoom = 4 for Z=21.
%% **********************************************************************
load Hb_ram;
load Hb_1112;

B = setup.Blength;
R = setup.Rate;
flag = setup.flag;

switch R
    case 7/8, Hb = Hb78;
    case 3/4, Hb = Hb34;
    case 1/2
        if flag
            switch B
                case 8064, Hb  = Hb_8k_1112;
                case 4032, Hb  = Hb_4k_1112;
                case 2016, Hb  = Hb_2k_1112;
                otherwise, error('Block length not supported')
            end
        else
            Hb = Hb12;
        end
    otherwise, error('Rate not supported')
end

switch B
    case  8064, Z = 84; zz = 1;
    case  4032, Z = 42; zz = 2;
    case  2016, Z = 21; zz = 4;
end

if flag
    if R==3/4 || R == 7/8
        for n=1:size(Hb,1)
            for m=1:size(Hb,2)
                if Hb(n,m)~=-1
                    Hb(n,m)=floor(Hb(n,m)/zz);
                end
            end
        end
    end
else
    
    for n=1:size(Hb,1)
        for m=1:size(Hb,2)
            if Hb(n,m)~=-1
                Hb(n,m)=floor(Hb(n,m)/zz);
            end
        end
    end
    
end
[M,N] = size(Hb);
Hm = logical(zeros(M*Z,N*Z));

for m=1:M
    for n=1:N
        if Hb(m,n)~=-1
           Hm(Z*(m-1)+1:m*Z,Z*(n-1)+1:n*Z) = circshift(logical(eye(Z)),Hb(m,n))';
        end
    end
end
Hx = sparse(Hm);
Hobj.Hbin = Hm;
Hobj.Hsprx = Hx;
Hobj.Rate = R;
Hobj.Blength = B;
end