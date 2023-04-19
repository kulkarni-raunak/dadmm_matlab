
function Hhat = symposdef(H)

% In order to check if a matrix is symmetric positive definite matrix if
% not then convert it into a symmetric positive definite matrix


% FIrst test for a square matrix A
[r,c] = size(H);
if r ~= c
  error('A must be a square matrix.')
elseif (r == 1) && (H <= 0)   
  % A was scalar and non-positive, so just return eps
  Hhat = eps;
  return
end
% symmetrize H into B
B = (H + H')/2;
% Compute the symmetric polar factor of B. Call it H.
% Clearly H is itself SPD.
[U,Sigma,V] = svd(B);
H = V*Sigma*V';
% get Hhat
Hhat = (B+H)/2;
% ensure symmetry
Hhat = (Hhat + Hhat')/2;
% test that Hhat is in fact PD. if it is not so, then tweak it just a bit.
p = 1;
k = 0;
while p ~= 0
  [R,p] = chol(Hhat);
  k = k + 1;
  if p ~= 0
    mineig = min(eig(Hhat));
    Hhat = Hhat + (-mineig*k.^2 + eps(mineig))*eye(size(H));
  end
end
