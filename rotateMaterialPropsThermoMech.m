function [Dcomp, Kacomp, Acomp, dDlcomp, dDmcomp, dDncomp, dDycomp, ...
          dKaycomp, dKalcomp, dKamcomp, dKancomp, dAycomp, dAlcomp, dAmcomp, dAncomp] = ...
          rotateMaterialPropsThermoMech(D0, Ka0, A0, l, m, n, dD0y, dKa0dy, dA0dy)
% rotateMaterialPropsThermoMech
% -------------------------------------------------------------------------
% Computes the 21 lower-triangular components of the rotated 3D stiffness 
% matrix, the 6 lower-triangular components of the rotated thermal 
% conductivity matrix, and the 6 components of the rotated thermal 
% expansion vector, along with their directional and content derivatives.
%
% This function fully vectorizes the calculations using compact lower-triangular 
% representations to optimize memory usage and performance, avoiding heavy 
% 3D tensor/page-wise matrix operations.
%
% Voigt Notation Mappings:
%   Dcomp  (Stiffness)          - 21 lower-triangular components (Voigt order)
%   Kacomp (Thermal Cond)       - 6 lower-triangular components (3x3 symmetric: [11, 21, 22, 31, 32, 33])
%   Acomp  (Thermal Exp Vector) - 6 Voigt vector components ([11, 22, 33, 23, 13, 12])
%
% Inputs:
%   D0, Ka0, A0 - Local material matrices (6x6xN, 3x3xN, 3x3xN)
%   l, m, n     - Direction cosines of fiber orientation (Nx1)
%   dD0y, dKa0dy, dA0dy - Derivatives w.r.t fiber volume fraction y (Optional)
%
% Outputs:
%   Dcomp, Kacomp, Acomp - Rotated matrices/vectors in vectorized format
%   dDlcomp, ..., dAncomp - Direct directional derivatives (w.r.t l, m, n)
%   dDycomp, dKaycomp, dAycomp - Content fraction derivatives (w.r.t y)
% -------------------------------------------------------------------------

% Ensure column vectors for direction cosines
l = l(:);
m = m(:);
n = n(:);

% =========================================================================
% 1. MECHANICAL STIFFNESS MATRIX ROTATION (21 components)
% =========================================================================
% Extract independent local stiffness components
C11_D = squeeze(D0(1,1,:)); C11_D = C11_D(:);
C22_D = squeeze(D0(2,2,:)); C22_D = C22_D(:);
C12_D = squeeze(D0(1,2,:)); C12_D = C12_D(:);
C23_D = squeeze(D0(2,3,:)); C23_D = C23_D(:);
C44_D = squeeze(D0(5,5,:)); C44_D = C44_D(:);

% Convert to invariant stiffness coefficients
a1 = C23_D;
a2 = (C22_D - C23_D) / 2;
a3 = C12_D - C23_D;
a4 = C44_D - (C22_D - C23_D) / 2;
a5 = C11_D + C22_D - 2*C12_D - 4*C44_D;

% Build rotated 21 lower-triangular mechanical stiffness components
Dcomp = fillStiffnessLowerTri21(a1, a2, a3, a4, a5, l, m, n);

% Calculate stiffness directional derivatives
[dDlcomp, dDmcomp, dDncomp] = computeDirectionalDerivsLowerTri21(a3, a4, a5, l, m, n);

% Calculate stiffness derivative with respect to fiber content y
if nargin >= 7 && ~isempty(dD0y)
    dC11y = squeeze(dD0y(1,1,:)); dC11y = dC11y(:);
    dC22y = squeeze(dD0y(2,2,:)); dC22y = dC22y(:);
    dC12y = squeeze(dD0y(1,2,:)); dC12y = dC12y(:);
    dC23y = squeeze(dD0y(2,3,:)); dC23y = dC23y(:);
    dC44y = squeeze(dD0y(5,5,:)); dC44y = dC44y(:);

    da1y = dC23y;
    da2y = (dC22y - dC23y) / 2;
    da3y = dC12y - dC23y;
    da4y = dC44y - (dC22y - dC23y) / 2;
    da5y = dC11y + dC22y - 2*dC12y - 4*dC44y;

    dDycomp = fillStiffnessLowerTri21(da1y, da2y, da3y, da4y, da5y, l, m, n);
else
    dDycomp = [];
end

% =========================================================================
% 2. THERMAL CONDUCTIVITY MATRIX ROTATION (6 components)
% =========================================================================
% Extract independent local thermal conductivity components
C11_K = squeeze(Ka0(1,1,:)); C11_K = C11_K(:);
C22_K = squeeze(Ka0(2,2,:)); C22_K = C22_K(:);

% Difference coefficient
g = C11_K - C22_K;

% Build rotated 6 lower-triangular thermal conductivity components
Kacomp = fillThermalCondLowerTri6(C22_K, g, l, m, n);

% Calculate thermal conductivity directional derivatives
[dKalcomp, dKamcomp, dKancomp] = computeThermalCondDerivsLowerTri6(g, l, m, n);

% Calculate thermal conductivity derivative with respect to fiber content y
if nargin >= 8 && ~isempty(dKa0dy)
    dC11Ky = squeeze(dKa0dy(1,1,:)); dC11Ky = dC11Ky(:);
    dC22Ky = squeeze(dKa0dy(2,2,:)); dC22Ky = dC22Ky(:);
    gy = dC11Ky - dC22Ky;
    dKaycomp = fillThermalCondLowerTri6(dC22Ky, gy, l, m, n);
else
    dKaycomp = [];
end

% =========================================================================
% 3. THERMAL EXPANSION VECTOR ROTATION (6 components)
% =========================================================================
% Extract independent local thermal expansion components
C11_A = squeeze(A0(1,1,:)); C11_A = C11_A(:);
C22_A = squeeze(A0(2,2,:)); C22_A = C22_A(:);

% Difference coefficient
h = C11_A - C22_A;

% Build rotated 6 thermal expansion components (Voigt Vector format)
Acomp = fillThermalExpLowerTri6(C22_A, h, l, m, n);

% Calculate thermal expansion directional derivatives
[dAlcomp, dAmcomp, dAncomp] = computeThermalExpDerivsLowerTri6(h, l, m, n);

% Calculate thermal expansion derivative with respect to fiber content y
if nargin >= 9 && ~isempty(dA0dy)
    dC11Ay = squeeze(dA0dy(1,1,:)); dC11Ay = dC11Ay(:);
    dC22Ay = squeeze(dA0dy(2,2,:)); dC22Ay = dC22Ay(:);
    hy = dC11Ay - dC22Ay;
    dAycomp = fillThermalExpLowerTri6(dC22Ay, hy, l, m, n);
else
    dAycomp = [];
end

end