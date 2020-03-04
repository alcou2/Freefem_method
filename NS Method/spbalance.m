function [T, B] = spbalance(A)
    % SPBALANCE balancing for sparse matrices.
    %   SPBALANCE attempts to bring the row and column norms of a square
    %   matrix closer together with a similarity transformation,
    %   consisting of a permutation and a diagonal scaling.  The scaling
    %   uses powers of two to avoid rounding errors. Balancing a matrix
    %   may improve the accuracy of computed eigenvalues.
    %
    %   B = spbalance(A) returns the balanced matrix B.
    %
    %   [T, B] = spbalance(A) returns the transformation matrix T and
    %   the balanced matrix B such that B = T\A*T.
    %
    %   To use SPBALANCE, the MEX-file "mexspbalance.c" must be
    %   compiled first. This can be done by executing the command "mex
    %   -largeArrayDims mexspbalance.c".
    %
    %   See also BALANCE.
    %
    %   Copyright Ian Zwaan <i.n.zwaan@tue.nl> 2013 - 2016.  Distributed
    %   under the Boost Software License, Version 1.0; see
    %   mexspbalance.c or copy at http://www.boost.org/LICENSE_1_0.txt.

    narginchk(1, 1);
    nargoutchk(1, 2);

    assert(issparse(A));
    assert(issquare(A));
    assert(isa(A, 'double'));
    assert(all(isfinite(A(find(A)))));

    if 0 == nnz(A)
        B = A;
        T = speye(size(A));
    else
        [T, B] = mexspbalance(A);
    end

    if 1 == nargout
        T = B;
    end
end

function b = issquare(A)
    b = size(A, 1) == size(A, 2);
end
