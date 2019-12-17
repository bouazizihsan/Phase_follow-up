function OUT=create_constellation1(m,MAPPING,MOD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% m: number of bits per symbol
%%% MAPPING: % 'Gray', 'Ungerboeck'
%%% OUT.X - constellation
%%% OUT.Lbin - binary mapping
%%% OUT.pntk0 - pointers to the symbols labeled with 0
%%% OUT.pntk1 - pointers to the symbols labeled with 1
%%% OUT.Xsort - constellation sorted in ascending order of the decimal
%%%             labels, so that OUT.Xsort( d ) is the symbol with the decimal label d
%%% author: LS, Fev. 2017

M=2^m;
switch MOD
    case 'QAM'
        % create QAM constellation
        if 2*floor(m/2)==m  % m is even
            pam=[-(sqrt(M)-1):2:sqrt(M)-1];     	% PAM constellation
            [A,B]=meshgrid(pam,pam(end:-1:1));   	% QAM
            X=A+1i*B;X=X(:);
            X=X/sqrt(mean(abs(X).^2));
            dx=real(X(2))-real(X(1));
            
            % define the labeling
            switch MAPPING
                case 'Gray', %% Gray
                    %%%% Create the Gray labelling for M-QAM
                    LPAM 	= Get_Labeling(m/2,'BRGC'); % Get the labeling (for PAM)
                    L       = bi2de(LPAM,'left-msb');
                    [A,B]  	= meshgrid(L,L);                            % Create the labeling for QAM
                    CAC     = bi2de([de2bi(A(:),m/2,'left-msb'),de2bi(B(:),m/2,'left-msb')],'left-msb');
                    Lbin    = de2bi(CAC,m,'left-msb');
                    %%% sort : the labels are in ascending order
                    [zz,fs]= sort(CAC(:));
                    Xsort       = X(fs);  %% reshuffle the constellation correspondingly
                    Lbin=Lbin(fs,:);
                    for kk=1:m                      % Create vectors with pointers to subconstellations
                        pntk1(:,kk)=find(Lbin(:,kk)==1);
                        pntk0(:,kk)=find(Lbin(:,kk)==0);
                    end
                case  'Ungerboeck', %%Ungerboeck with 16x M/16 partition, for M>=64
                    if M<64, error('Ungerboeck mapping applicable only for 64QAM and higher'); end
                    %%%% First Create the Gray labelling for 16-QAM
                    LPAM 	= Get_Labeling(2,'BRGC'); % Get the labeling (for 16 PAM)
                    L       = bi2de(LPAM,'left-msb');
                    [A,B]  	= meshgrid(L,L);                            % Create the labeling for QAM
                    CAC16   = reshape(bi2de([de2bi(A(:),2,'left-msb'),de2bi(B(:),2,'left-msb')],'left-msb'),4,4);
                    % repeat the Gray mapping to form the Ungerboeck mapping
                    CAC     = repmat(CAC16,sqrt(M)/4,sqrt(M)/4);
                    % --Create the Gray mapping for (m-4) uncoded bits
                    LPAM 	= Get_Labeling((m-4)/2,'BRGC');
                    L       = bi2de(LPAM,'left-msb');
                    CAC_u   = zeros(sqrt(M),sqrt(M));
                    tt=(1:4);t_off=4;
                    for ik=0:sqrt(M)/4 -1
                        for il=0:sqrt(M)/4 -1
                            CAC_u( tt+ik*t_off , tt+il*t_off ) = L(ik+1)*2^(m/2-2) +L(il+1);
                        end
                    end
                    % --
                    % Combine repeated Gray mapping for Ungerboeck labels with the
                    % label of uncoded bits
                    CAC     = CAC*( 2^(m-4) ) + CAC_u;
                    Lbin    = de2bi(CAC(:),m,'left-msb');
                    %%% sort : the labels are in ascending order
                    [zz,fs]= sort(CAC(:)); % A=CAC(:); zz= A(fs)
                    Xsort       = X(fs);  %% reshuffle the constellation correspondingly
                    Lbin=Lbin(fs,:);
                    
                    for kk=1:m    % Create vectors with pointers to subconstellations
                        pntk1(:,kk)=find(Lbin(:,kk)==1);
                        pntk0(:,kk)=find(Lbin(:,kk)==0);
                    end
                case 'Bin',
                    X=qammod(0:M-1,M,0,'bin');
                    power=sum(abs(X).^2)./length(X);
                    X=X./sqrt(power);
                    X=X.';
                    Xsort=X;
                    [pntk1,pntk0]=testbit(M);
                    
                    
            end
            OUT.X = X;
            
            OUT.pntk0 = pntk0;
            OUT.pntk1 = pntk1;
            OUT.Xsort = Xsort;
            OUT.MAPPING = MAPPING;
            OUT.phi_min = atan(1/pam(end-1))-atan(1/pam(end));
            
        else %% m is odd, Cross-shaped constellations ( m > 3 )
            %             M2=M/2;
            %             pam=[-(sqrt(M2)-1):2:sqrt(M2)-1];     	% PAM constellation
            %             [A,B]=meshgrid(pam,pam(end:-1:1));   	% inner QAM
            %             X=A+1i*B;
            %             K=sqrt(M2)/4;  %% crosses arm size
            %             Cnorth=X(1:K,:)+1i*K*2 ;
            %             Csouth=X(end-K+1:end,:)-1i*K*2;
            %             Cwest=X(:,1:K)-K*2;
            %             Ceast=X(:,end-K+1:end)+K*2;
            %             X=[X(:);Cnorth(:);Csouth(:);Cwest(:);Ceast(:)];
            %             X=X/sqrt(mean(abs(X).^2));
            X=qammod(0:M-1,M,0,'gray');
            power=sum(abs(X).^2)./length(X);
            X=X./sqrt(power);
            OUT.X=X.';
            OUT.Xsort=X.';
            [pntk1,pntk0]=testbit(M);
            OUT.pntk0 = pntk0;
            OUT.pntk1 = pntk1;
            OUT.MAPPING = MAPPING;
        end
    case 'PSK'
        X=pskmod(0:M-1,M);
        power=sum(abs(X).^2)./length(X);
        X=X./sqrt(power);
        OUT.X=X.';
        OUT.Xsort=X.';
        [pntk1,pntk0]=testbit(M);
        OUT.pntk0 = pntk0;
        OUT.pntk1 = pntk1;
        OUT.MAPPING = MAPPING;
    otherwise
        error('unknown Modulation Type');
end
end
function L=Get_Labeling(m,type)
% L=Get_Labeling(m,type) returns the decimal representation of some
% well-known binary labelings. L is a length-M column vector and m is the
% number of bits per symbol, i.e., m=log2(M);
% The supported labelings are 'BRGC', 'NBC', and 'AGC'
%
% Code originally written by Fredrik Br?nnstr?m.
% Modified by Alex Alvarado, June 2014
% Modified by LS, Fev 2017

if m==1
    L=[0 1]';
else
    M=2^m;
    switch type
        case 'BRGC'
            L=zeros(M,m);
            L(1:M/2,2:m)=Get_Labeling(m-1,type);
            L(M/2+1:M,2:m)=flipud(L(1:M/2,2:m));
            L(M/2+1:M,1)=1;
        case 'NBC'
            L=fliplr(de2bi(0:M-1));
        case 'AGC'  % This is the AGC for PAM
            L=zeros(M,m);
            L(1:M/2,2:m)=Get_Labeling(m-1,type);
            L(M/2+1:M,2:m)=flipud(L(1:M/2,2:m));
            L(2:2:M,1)=1;
            L(M/2+1:M,:)=not(L(M/2+1:M,:));
        case 'UNG' %% Ungerboeck labelling
            
        otherwise
            error(sprintf('Only ''BRGC'', ''NBC'', and ''AGC'' are supported and here type=''%s''',type))
    end
end
end

function [outbit1,outbit0]=testbit(n)
m=log2(n);
outbit1=zeros(n/2,m);
outbit0=zeros(n/2,m);
for k=1:m
    b=bitget(0:n-1,m-k+1);% bitget(nombre decimal, bit), nombre decimal se transforme automatiquement en binaire, right-msb donc (m-k+1)
    [row,col] = find(b);
    outbit1(:,k)=col;
    [row,col] = find(b==0);
    outbit0(:,k)=col;
end
end


