function e = cartesianbasis2d(e1, e2)

% CARTESIANBASIS2D Construct two-dimensional cartesian basis vectors.
%   CARTESIANBASIS2D returns a column matrix containing two unit vectors
%   which form an orthonormal basis. The basis thus defined may
%   subsequently be used to define vectors and tensors.
%
%   CARTESIANBASIS2D(s1,s2) names the two basis vectors s1 and s2
%   respectively. If (one or both of) these arguments are absent, the
%   default names e(1), e(2) are used.
%
%   Example
%      e = cartesianbasis2d('ex', 'ey');
%      ex = e(1);
%      ey = e(2);
%      a = ex + 2*ey
%      B = ex*ex + ex*ey - ey*ex
%
%   See also CARTESIANBASIS1D, CARTESIANBASIS3D.

% complete arguments
if (nargin < 1) || isempty(e1)
    e1 = 'e(1)';
end
if (nargin < 2) || isempty(e2)
    e2 = 'e(2)';
end

% construct structure
e.basis = str2mat(e1, e2);
e.order = 1;
e.size  = [2 1];
e.components = [1 0 0 1]';

% form tensor object
e = tensor(e);

end