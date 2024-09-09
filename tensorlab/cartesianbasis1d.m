function e = cartesianbasis1d(e1)

% CARTESIANBASIS1D Construct one-dimensional cartesian basis vector.
%   CARTESIANBASIS1D returns a unit vector which can be used as a
%   one-dimensional basis. This basis may subsequently be used to define
%   one-dimensional vectors and tensors.
%
%   CARTESIANBASIS1D(s) names the basis vector s. If this argument is
%   absent, the default name e(1) is used.
%
%   Example
%      ex = cartesianbasis2d('ex');
%      a = ex
%      B = ex*ex
%
%   See also CARTESIANBASIS2D, CARTESIANBASIS3D.

if (nargin == 0) || isempty(e1)
    e1 = 'e(1)';
end
e.basis = e1;
e.order = 1;
e.size  = [1 1];
e.components = 1;

e = tensor(e);
