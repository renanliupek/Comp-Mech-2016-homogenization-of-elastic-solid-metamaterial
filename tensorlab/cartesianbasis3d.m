function e = cartesianbasis3d(e1, e2, e3)

% CARTESIANBASIS3D Construct three-dimensional cartesian basis vectors.
%   CARTESIANBASIS3D returns a column matrix containing three unit vectors
%   which form a right-handed, orthonormal basis. The basis thus defined
%   may subsequently be used to define vectors and tensors.
%
%   CARTESIANBASIS3D(s1,s2,s3) names the three basis vectors s1, s2 and s3
%   respectively. If (one or several of) these arguments are absent, the
%   default names e(1), e(2) and e(3) are used.
%
%   Example
%      e = cartesianbasis3d('ex', 'ey', 'ez');
%      ex = e(1);
%      ey = e(2);
%      ez = e(3);
%      a = ex + 2*ey + 3*ez
%      B = ex*ex + ex*ey - ey*ex + 2*ez*ez
%
%   See also CARTESIANBASIS1D, CARTESIANBASIS2D.

% complete arguments
if (nargin < 1) || isempty(e1)
    e1 = 'e(1)';
end
if (nargin < 2) || isempty(e2)
    e2 = 'e(2)';
end
if (nargin < 3) || isempty(e3)
    e3 = 'e(3)';
end

% construct structure
e.basis = str2mat(e1, e2, e3);
e.order = 1;
e.size  = [3 1];
e.components = [1 0 0 0 1 0 0 0 1]';

% form tensor object
e = tensor(e);

end