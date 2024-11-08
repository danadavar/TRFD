function [A, H, nf] = Jac_approx(x, fvec, nf, tau, m, n, hfun, Ffun)

H = [];
I = eye(n);
A = zeros(m,n);

for j = 1:n
  z = x+tau*I(:,j);
  fvecz = Ffun(z);
  nf = nf+1;
  H = [H;hfun(fvecz)];
  A(:,j) = (fvecz-fvec)/tau;
end
