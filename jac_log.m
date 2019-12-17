function z=jac_log(x,y)

%%%%%%%%%%%
% jacobian lagorithm to calculate log( sum( exp(x) ) )
if nargin==2  %%% put both arguments as separate columns
    z=jac_log( [x(:) y(:)] );
    return
end

%%%% if x is a matrix, jacobian logarithm is applied to each row 

%%% the new non-recurive version
K=size(x,2);
z=max(x,[],2);
z=z+log(sum(exp(x-repmat(z,1,K)),2));

%%%% this is the old recursive version
% switch size(x,2)
%    case 1, z=x;
%    case 2, z=max(x,[],2)+log( 1+exp(-abs( x(:,2)-x(:,1) ) ) ); 
%    otherwise, %>2
%     y=jac_log( x(:,1:end-1) ); %% recursive formulation
%     z=jac_log([y x(:,end)]);
%     return
% end

