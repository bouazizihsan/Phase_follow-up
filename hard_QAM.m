function [xhat, ii] = hard_QAM( y , X )

%   hard decisions on y 
%   X is a QAM constellation in listed column-wise starting from
%      the upper-left corner

M = length(X);
K = sqrt(M);
D = abs( X(2)-X(1) );

kx = min( max( floor( real(y)./D ) + K/2, 0) , K-1) ; % column index
ky = min( max( floor( imag(y)./D ) + K/2, 0), K-1); % row index

ii = (kx+1)*K - ky;
xhat = X( ii );

