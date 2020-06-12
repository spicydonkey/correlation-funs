function dtheta = pair_diff_angle(X1,X2)
% pair_diff_angle finds all pairwise difference angles between Cart vectors
% in X1 and X2.
%   
% dtheta = pair_diff_angle(X1,X2)   pairs X1 and X2
% dtheta = pair_diff_angle(X1)      pairs with itself
% 
% X1: Nx3 array of 3D cart vectors
% X2: Mx3 array of 3D cart vectors
%
% dtheta: relative angle between vector in X1 and X2
% 
% NOTE: dtheta for self is set to NaN if self-pairing occurs;
%
% DKS 2020

b_self = false;
if nargin<2
    b_self = true;
    X2 = X1;
end

[x1,y1,z1] = columns(X1);
[x2,y2,z2] = columns(X2);

dtheta = diffAngleXYZ(x1,y1,z1,x2',y2',z2');

if b_self
    dtheta(boolean(eye(size(dtheta)))) = nan;
end

end