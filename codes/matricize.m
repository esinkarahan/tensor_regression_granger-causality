function Tm = matricize(T,varargin) 
%   Tm = TENMAT(T, RDIMS, CDIMS) creates a matrix representation of
%   tensor T.  The dimensions specified in RDIMS map to the rows of
%   the matrix, and the dimensions specified in CDIMS map to the
%   columns, in the order given.
% Modified tenmat function of Tensor Toolbox of Kolda et. al.
% MATLAB Tensor Toolbox.
% Copyright 2010, Sandia Corporation. 

    tsize = size(T);
    n = ndims(T);
    if nargin==2
        rdims=varargin{1};
        tmp = true(1,n);
        tmp(rdims) = false;
        cdims = find(tmp);
    else
        cdims = varargin{1};
        tmp = true(1,n); 
        tmp(cdims) = false; 
        rdims = find(tmp);
    end

    data = reshape(permute(T,[rdims cdims]), prod(tsize(rdims)), prod(tsize(cdims)));
    Tm.tsize = tsize;
    Tm.rindices = rdims;
    Tm.cindices = cdims;
    Tm.data = data;
end


