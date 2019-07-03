function [opts,alphaL1,alphaL2,alphaorth,keepsign,LtL,L] = setparametersTF_Shorthals (opts,N,szx)
% Set parameters for the penalized Tensor Factorization
% INPUTS
% opts   : structure to be set
% N        : order of the tensor
% szx      : size of the tensor
% 
% OUTPUTS
% 
% opts structure has cell fields specifying the penalization parameters for
% each factor of the decomposition, these are
% opts.alphaorth for the orthogonality parameter
% opts.alphaL1 for the sparsity parameter
% opts.alphaL2 for the smoothness parameter
% Additionaly the algorithmOpts field is a cell array contatining the
% laplacian smoother matrices for each factor
%
% Version 1 - May 2015 
%

if ~isfield(opts,'tol'),         opts.tol = 1e-4;           end
if ~isfield(opts,'maxiters'),    opts.maxiters = 100;       end
if ~isfield(opts,'verbose'),     opts.verbose = 1;          end
if ~isfield(opts,'init'),        opts.init = 'random';      end
if ~isfield(opts,'const'),       opts.const = zeros(N,1);   end
if ~isfield(opts,'makezero'),    opts.makeZero = eps;      end % zero the values which are near to 0

% Find the factor that will keep the sign for real tensor factorization
% The factor which is real
keepsign = find(~opts.nonnegative);
try      keepsign = keepsign(1);
catch,   keepsign = 1; end

if ~isfield(opts,'algorithmOpts'), 
    for n = 1:N, opts.algorithmOpts{n} = struct; end        
end

% Algorithm options
L   = cell(N,1);
LtL = cell(N,1);
alphaL1   = cell(N,1);
alphaL2   = cell(N,1);
alphaorth = cell(N,1);

for n = 1:N
    if opts.const(n)== 1 % Orthogonality + Smooth Lasso
        
        alphaL1{n} = setAlpha(opts,'alphaL1',n);
        alphaL2{n} = setAlpha(opts,'alphaL2',n);
        alphaorth{n} = setAlpha(opts,'alphaorth',n);
        
        if ~isfield(opts.algorithmOpts{n},'laplace'), opts.algorithmOpts{n}.laplace = 1; end
              
        L = findL(L,opts,n,szx);
        LtL{n} = L{n}'*L{n};
        
    elseif opts.const(n)== 2 % Smooth Lasso
        
        alphaL1{n} = setAlpha(opts,'alphaL1',n);
        alphaL2{n} = setAlpha(opts,'alphaL2',n);
        
        if ~isfield(opts.algorithmOpts{n},'laplace'), opts.algorithmOpts{n}.laplace = 1; end

        L   = findL(L,opts,n,szx);
        LtL{n} = L{n}'*L{n};

    end
end

   

end


function L = findL(L,opts,n,szx)
    if isfield(opts.algorithmOpts{n},'L'),
        L{n} = opts.algorithmOpts{n}.L;
    else
        if isfield(opts.algorithmOpts{n},'laplace') % add laplacian
            if opts.algorithmOpts{n}.laplace == 1
                % 1D Laplacian
                e = ones(szx(n), 1);
                L{n} = spdiags([-e e], 0:1, szx(n), szx(n));
            else
                % 2D Laplacian
                e = ones(szx(n), 1);
                L{n} = spdiags([e -2*e e], -1:1, szx(n), szx(n));
            end
        else
            L{n} = speye(szx(n));
        end
    end
end
%--------------------------------------------------------------------------

function alpha = setAlpha(opts,alphaType,n)
if ~isfield(opts,alphaType),
    alpha = 1;
else
    tmp = getfield(opts,alphaType);
    alpha=tmp{n};
end
end

