function [x, f_min, nf, stop, H] = TRFD_composite(x0, m, Ffun, h, nfmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%
% This function attempts to solve composite problems of the form
%
% min f(x) = h(F(x)), F: R^(n) -> R^(m), h: R^(m) -> R,
% where h is either the L1-norm or the Max-function,
%
% by applying the method proposed in:
%
% Davar, D., Grapiglia, G. N.: TRFD: A derivative-free trust-region method 
% based on finite differences for composite nonsmooth optimization
%
% Input:
%  
%   x0 [n x 1] = initial point (a column vector)
%   m = number of components of the inner function F
%   Ffun = function provided by the user that computes F
%   h [1 x 1] = a scalar defining the outer function h:
%               1 for L1 problems,
%               2 for Minimax problems
%   nf_max = maximum number of function evaluations allowed
%
% Output:
%
%   f_min = lowest function value found within nfmax function evaluations
%   x     = vector yielding the lowest function value: f_min = f(x)
%   nf    = number of function evaluations used to obtain x 
%   stop  = an integer identifying the reason the algorithm stopped:
%           
%        stop = 1 : Number of function evaluations >= nf_max
%        stop = 2 : Trust-region radius less than Delta_tol (default value = 1e-13)
%        stop = 3 : Approximate criticality measure less than eta_tol (default value = 1e-13)
%        stop = -1 : Error in the execution of linprog
%        
%    H = a vector such that H(i) is the smallest value
%        of the objective f obtained by the algorithm after 
%        i function evaluations.
%
% Functions called: Ffun, Jac_approx, inner_solver_L1, inner_solver_max
%                   
%
% Authors information:
%
% Dânâ Davar, Geovani Nunes Grapiglia
% Université catholique de Louvain, Belgium
% Department of Mathematical Engineering 
% dana.davar@uclouvain.be, geovani.grapiglia@uclouvain.be
%
% November 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long

%%%%%%%%%%%%%%%%%%%  OUTER FUNCTION AND NORM-CONSTANTS %%%%%%%%%%%%%%%%%%%%

n = height(x0);                  % number of variables

if h == 1
    hfun = @(z) norm(z,1);
    L_h = 1;                     % Lipschitz constant of the outer function
    cm = sqrt(m);                % ||z||_1 <= cm ||z||_2
    cn = 1;                      % ||x||_2 <= cn ||x||_1
end

if h == 2
    hfun = @(z) max(z);
    L_h = 1;                    % Lipschitz constant of the outer function
    
    if sqrt(m) >= n
        cm = 1;                 % ||z||_inf <= cm ||z||_2                         
        cn = sqrt(n);           % ||x||_2 <= cn ||x||_inf
    else
        cm = sqrt(m);           % ||z||_1 <= cm ||z||_2                  
        cn = 1;                 % ||x||_2 <= cn ||x||_1
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%  ALGORITHMIC PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%

ep = 10^(-15);                              % tolerance epsilon
sigma = ep/(L_h*cm*cn*sqrt(n)*sqrt(eps));   % sigma
tau = ep/(L_h*sigma*cm*cn*sqrt(n));         % finite-difference step-size
alpha = 0.15;                               % minimum ratio to label the iteration as successful
Delta = max([1 tau*sqrt(n)]);               % initial trust-region radius
Delta_max = 1000;                           % upper-bound on the trust-region radius
Delta_tol = 1e-13;                          % tolerance on the trust-region radius
eta_tol = 1e-13;                            % tolerance on the approximate criticality measure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  STEP 0  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = x0;                                     % initial point x0
fvec = Ffun(x);                             % F(x0)
f = hfun(fvec);                             % h(F(x0))
nf = 1;                                     % number of function evaluations
H0 = [f];                                   % history
k = 0;                                      % iteration index
options = optimset('Display','none');       % linprog option to avoid displaying
stop = 0;                               

Step = 1;

while k >= 0 && nf < nfmax && Delta > Delta_tol

    switch Step

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  STEP 1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        case 1
    
            [A, H1, nf] = Jac_approx(x, fvec, nf, tau, m, n, hfun, Ffun);
            H0 = [H0;H1];
          
            if h == 1
                [d, exit] = inner_solver_L1(A, m, n, fvec, Delta_max, options);
            end

            if h == 2
                [d, exit] = inner_solver_max(A, m, n, fvec, Delta_max, options);
            end
          
            if exit == -1
               stop = -1;
               break
            end
          
            eta = (f-hfun(fvec+A*d))/Delta_max;
    
            if eta <= eta_tol
               stop = 3;
               break
            end
    
            Step = 2;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  STEP 2  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        case 2
    
            if eta >= 0.5*ep
                Step = 3;
            else
                tau = tau/2;
                k = k+1;
                Step = 1;
            end
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  STEP 3  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        case 3
    
            if h == 1
                [d, exit] = inner_solver_L1(A, m, n, fvec, Delta, options);
            end

            if h == 2
                [d, exit] = inner_solver_max(A, m, n, fvec, Delta, options);
            end

            if exit == -1
               stop = -1;
               break
            end
      
            Step = 4;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  STEP 4  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        case 4
    
            x1 = x+d;
            fvec1 = Ffun(x1);
            f1 = hfun(fvec1);
            nf = nf+1;
            H0 = [H0;f1];

            Ared = f - f1;
            Pred = f - hfun(fvec+A*d); 
        
            rho = Ared/Pred;
    
            if rho >= alpha
                Delta = min([2*Delta Delta_max]);
                x = x1;
                fvec = fvec1;
                f = f1;
                k = k+1;
                Step = 1;
            else
                Step = 5;
            end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  STEP 5  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        case 5
    
            Delta = Delta/2;
    
            if tau*sqrt(n) <= Delta
                k = k+1;
                Step = 3;
            else
                tau = tau/2;
                k = k+1;
                Step = 1;
            end

%%%%%%%%%%%%%%%%%%%%%%%%%%  SWITCH STEP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%  TRFD ENDS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

if nf >= nfmax
    stop = 1;
end

if Delta <= Delta_tol
    stop = 2;
end

%%%%%%%%%%%%%%%%%%%%%  BUILD MONOTONE HISTORY H  %%%%%%%%%%%%%%%%%%%%%%%%%%

if height(H0) == 1
    H0 = [H0;H0];
end

H = zeros(nfmax,1);
H(1) = H0(1);

if nfmax <= height(H0)
    for i = 2:nfmax
        if H0(i) > H(i-1)
            H(i) = H(i-1);
        else
            H(i) = H0(i);
        end
    end
else
    for i = 2:height(H0)
        if H0(i) > H(i-1)
            H(i) = H(i-1);
        else
            H(i) = H0(i);
        end
    end

    for j = height(H0)+1:nfmax
        H(j) = H(height(H0));
    end
end

f_min = H(end);
plot(H(1:nf));
xlabel('Number of evaluations'); 
ylabel('Function value');
