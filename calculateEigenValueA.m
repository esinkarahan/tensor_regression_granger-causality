function opts = calculateEigenValueA(A,szy,Ny,opts)
% Calculates the Eigenvalue Decomposition of the N dimensional tensor A,
% asuming it is used as a design matrix in a regression problem. 
% PROPACK package from http://sun.stanford.edu/~rmunk/PROPACK/ is used for
% fast SVD calculation.
% INPUTS
% A    : N dimensional tensor to be decomposed
% szy  : size of the observation tensor
% Ny   : Number of dimensions of the observation tensor
% opts : several parameters for the algorithm, check
%        tensorRegressionCPpenSL.m for more info
% OUTPUTS
% opts: the resulting eigenvectors and eigenvalues are stored in the fields 
%        A.V and A.ev respectively
%
% Version 1 - May 2015 

% This is required to use Propack functions 
global AA
global Afunc;
% Find the fastest way for matrix multiplication in Matlab
AA = sprand(10,10,0.1);
forwardType = findBestMultiply(AA,.2);
Afunc       = sprintf('AforwardTen_%d',forwardType);

cmoda = opts.cmoda;
cmodx = opts.cmodx;

sza = size(A); Na = ndims(A);

Nx = length(cmodx) + (Ny-(Na-length(cmoda)));
szx(cmodx) = sza(cmoda);
szx((cmodx(end)+1):Nx) = setdiff(szy,sza(setdiff(1:Na,cmoda)));

cmodat = length(cmoda)+1:Na; % contraction modes for A'A
cmody  = 1:length(cmodat);   % contraction modes for A'Y

At = permute(A,[cmoda setdiff(1:Na,cmoda)]);

if prod(sza(cmody))>=prod(sza(cmoda))
    AtA = contractTensor(At,A,cmodat,cmody,'m');
else
    AtA = contractTensor(A,At,cmoda,cmodx,'m');
end

AA = AtA;
[V,ev]  = laneig(Afunc,size(AA,1),size(AA,1),'LM'); 

% To use Matlab builtin function eig replace laneig with this line
% [V,ev] = eig(AtA);
ev = diag(ev);
[ev,ind] = sort(ev,'descend');
V  = V(:,ind);
ev(ev<1e-14*ev(1))=0;
r  = sum(ev>0);
ev = ev(1:r);
V  = V(:,1:r);

opts.A.V = V;
opts.A.ev= ev;


end
