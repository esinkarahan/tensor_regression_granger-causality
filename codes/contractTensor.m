function C = contractTensor(a,b,adims,bdims,varargin)
% tensor tensor contraction
% C = TTT(a,b,adims,bdims) computes the contracted product of
%   tensors a and b in the dimensions specified by the row vectors
%   adims and bdims.  The sizes of the dimensions specified by adims
%   and bdims must match; that is, size(a,adims) must equal
%   size(b,bdims).
% opts: 't' result is returned as tensor  
%       'm' result is returned as matrix
% Modified ttt function of Tensor Toolbox of Kolda et. al.
% MATLAB Tensor Toolbox.
% Copyright 2010, Sandia Corporation. 

    if nargin >4
        opts = varargin{1};
    else
        opts = 't';
    end
    
    A = matricize(a,adims,'t');
    B = matricize(b,bdims);
    tsiz = [A.tsize(A.rindices) B.tsize(B.cindices)];
    rindices = 1:length(A.rindices);
    cindices = (1:length(B.cindices)) + length(A.rindices);
    C = A.data * B.data;

    if opts == 't'
        order = [rindices,cindices];
        data = reshape(C, [tsiz(order) 1 1]);
        if numel(order) >= 2
            C = ipermute(data,order);
        else
            C = data;
        end
    end
    
end

