function [y, yeq, grady, gradyeq] = quad_con (x, n, m, r, instance)

switch instance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INSTANCE 2  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 2

        y = x'*diag([ones(1,n), zeros(1,m)])*x + zeros(n+m,1)'*x - r*r;
        yeq = [];
            
        if nargout > 2
            grady(:,1) = 2*diag([ones(1,n), zeros(1,m)])*x + zeros(n+m,1);
        end
    
        gradyeq = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INSTANCE 5  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 5
        
        y = x'*diag([ones(1,n), 0])*x + zeros(n+1,1)'*x - r*r;
        yeq = [];
            
        if nargout > 2
            grady(:,1) = 2*diag([ones(1,n), 0])*x + zeros(n+1,1);
        end
    
        gradyeq = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INSTANCE 7  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 7

        y = x'*diag(ones(1,n))*x + zeros(n,1)'*x - r*r;
        yeq = [];
            
        if nargout > 2
            grady(:, 1) = 2*diag(ones(1,n))*x + zeros(n,1);
        end
    
        gradyeq = [];
end
