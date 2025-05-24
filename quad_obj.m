function [y, grady] = quad_obj (x, n, m, A, instance)

switch instance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INSTANCE 2  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 2
    
        y = [zeros(n,1); ones(m,1)]'*x;
        
        if nargout > 1
            grady = [zeros(n,1); ones(m,1)];
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INSTANCE 5  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 5
    
        y = [zeros(n,1); 1]'*x;
    
        if nargout > 1
            grady = [zeros(n,1); 1];
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INSTANCE 7  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 7
    
        y = A*x;
    
        if nargout > 1
            grady = A;
        end
end
