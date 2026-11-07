function [D0,dD0df,S0,K0,A0,dK0df,dA0df] = compositeMixingParams(km1,km2,am,Em,vm,kf1,kf2,af1,af2,Ef1,Ef2,vf12,vf23,Gf12,f)
% compositeMixingParams
% -------------------------------------------------------------------------
% Computes the local 3D constitutive stiffness matrix D0, thermal 
% conductivity matrix K0, and thermal expansion matrix A0 of a 
% transversely isotropic composite material, along with their derivatives
% with respect to the fiber volume fraction f.
%
% Inputs:
%   km1   - Matrix thermal conductivity in direction 1
%   km2   - Matrix thermal conductivity in directions 2 and 3
%   am    - Matrix coefficient of thermal expansion
%   Em    - Matrix Young's modulus
%   vm    - Matrix Poisson's ratio
%   kf1   - Fiber thermal conductivity in the longitudinal direction
%   kf2   - Fiber thermal conductivity in the transverse direction
%   af1   - Fiber longitudinal coefficient of thermal expansion
%   af2   - Fiber transverse coefficient of thermal expansion
%   Ef1   - Fiber Young's modulus in the longitudinal direction
%   Ef2   - Fiber Young's modulus in the transverse direction
%   vf12  - Fiber major Poisson's ratio
%   vf23  - Fiber transverse Poisson's ratio
%   Gf12  - Fiber longitudinal shear modulus
%   f     - Fiber volume fraction field, size: nelz x nely x nelx
%
% Outputs:
%   D0    - Local stiffness matrix, size: 6 x 6 x nElem
%   dD0df - Derivative of D0 with respect to f, size: 6 x 6 x nElem
%   S0    - Local compliance matrix, size: 6 x 6 x nElem
%   K0    - Local thermal conductivity matrix, size: 3 x 3 x nElem
%   A0    - Local thermal expansion matrix, size: 3 x 3 x nElem
%   dK0df - Derivative of K0 with respect to f, size: 3 x 3 x nElem
%   dA0df - Derivative of A0 with respect to f, size: 3 x 3 x nElem
%
% Notes:
%   The material is described in its local coordinate system. The first
%   material direction is the fiber direction. The mechanical, thermal 
%   conductivity, and thermal expansion properties are derived using 
%   established mixture rules for composite materials.
% -------------------------------------------------------------------------

% Number of elements in each direction and total number of elements
[nelz,nely,nelx] = size(f);
nElem = nelx*nely*nelz;

% Reshape f into page format so that each element corresponds to one page
% in the constitutive matrices.
f = reshape(f,1,1,nElem);

% Allocate mechanical compliance matrix and its derivative with respect to f
S0 = zeros(6,6,nElem);
dSdf = zeros(6,6,nElem);

% Allocate thermal conductivity matrix and its derivative with respect to f
K0 = zeros(3,3,nElem);
dK0df = zeros(3,3,nElem);

% Allocate thermal expansion matrix and its derivative with respect to f
A0 = zeros(3,3,nElem);
dA0df = zeros(3,3,nElem);

% =========================================================================
% MECHANICAL PROPERTIES
% =========================================================================

% Matrix shear modulus and fiber transverse shear modulus
Gm = Em/2/(1+vm);
Gf23 = Ef2/2/(1+vf23);

% Effective engineering constants based on the selected mixture rules
% E1:  longitudinal Young's modulus
% E2:  transverse Young's modulus
% G12: longitudinal shear modulus
% G23: transverse shear modulus
%
% The square-root mixture rules have an unbounded derivative at f = 0.
% Keep f = 0 physically admissible and regularize only a very small
% neighborhood of zero.  The C1 regularization below is exactly sqrt(f)
% for f >= eps_f, satisfies sqrt_f_reg(0) = 0, and matches both the value
% and first derivative of sqrt(f) at f = eps_f.
eps_f = 1e-6;
sqrt_f_reg = zeros(size(f));
dsqrt_f_reg = zeros(size(f));
mask_reg = (f < eps_f);
mask_orig = ~mask_reg;

sqrt_f_reg(mask_orig) = sqrt(f(mask_orig));
dsqrt_f_reg(mask_orig) = 0.5 ./ sqrt(f(mask_orig));

f_reg = f(mask_reg);
sqrt_f_reg(mask_reg) = 1.5*f_reg/sqrt(eps_f) ...
    - 0.5*f_reg.^2/(eps_f^(3/2));
dsqrt_f_reg(mask_reg) = 1.5/sqrt(eps_f) ...
    - f_reg/(eps_f^(3/2));

E1 = f*Ef1 + (1-f)*Em;
E2 = Em./(1-sqrt_f_reg*(1-Em/Ef2));
G12 = Gm./(1-sqrt_f_reg*(1-Gm/Gf12));
G23 = Gm./(1-sqrt_f_reg*(1-Gm/Gf23));

% Effective Poisson's ratios
v12 = f*vf12 + (1-f)*vm;
v23 = E2/2./G23 - 1;

% Assemble the local compliance matrix S0.
% Voigt notation order: [11, 22, 33, 23, 13, 12].
% The normal strain block is 3 x 3, and the shear terms are diagonal.
S0(1:3,1:3,:) = [ 1./E1,     -v12./E1,  -v12./E1; ...
                  -v12./E1,   1./E2,     -v23./E2; ...
                  -v12./E1,  -v23./E2,    1./E2 ];
S0(4,4,:) = 1./G23;
S0(5,5,:) = 1./G12;
S0(6,6,:) = 1./G12;

% Convert compliance matrix to stiffness matrix for every element/page
D0 = pageinv(S0);

% Derivatives of the effective engineering constants with respect to f
% dE1/df follows directly from the linear rule of mixtures.
dE1 = Ef1 - Em;

% dE2/df, dG12/df, and dG23/df use the exact derivative of the
% regularized square-root function above, preserving consistency between
% the forward material model and its analytical sensitivities.
dE2 = (Em*(1-Em/Ef2).*dsqrt_f_reg)./(1-sqrt_f_reg*(1-Em/Ef2)).^2;
dG12 = (Gm*(1-Gm/Gf12).*dsqrt_f_reg)./(1-sqrt_f_reg*(1-Gm/Gf12)).^2;
dG23 = (Gm*(1-Gm/Gf23).*dsqrt_f_reg)./(1-sqrt_f_reg*(1-Gm/Gf23)).^2;

% Derivatives of the effective Poisson's ratios with respect to f
dv12 = vf12 - vm;
dv23 = 0.5*(dE2.*G23 - E2.*dG23)./G23.^2;

% Assemble derivative of the compliance matrix with respect to f.
dSdf(1:3,1:3,:) = [ -1./(E1.^2).*dE1, ...
                    -1./E1.*dv12 + v12./E1.^2.*dE1, ...
                    -1./E1.*dv12 + v12./E1.^2.*dE1; ...
                    -1./E1.*dv12 + v12./E1.^2.*dE1, ...
                    -1./E2.^2.*dE2, ...
                    -1./E2.*dv23 + v23./E2.^2.*dE2; ...
                    -1./E1.*dv12 + v12./E1.^2.*dE1, ...
                    -1./E2.*dv23 + v23./E2.^2.*dE2, ...
                    -1./E2.^2.*dE2 ];

% Shear compliance derivatives
dSdf(4,4,:) = -1./G23.^2.*dG23;
dSdf(5,5,:) = -1./G12.^2.*dG12;
dSdf(6,6,:) = -1./G12.^2.*dG12;

% Derivative of stiffness matrix using D = inv(S):
%   dD/df = -D * (dS/df) * D
dD0df = pagemtimes(-pagemtimes(D0,dSdf),D0);

% =========================================================================
% THERMAL CONDUCTIVITY PROPERTIES
% =========================================================================

% Effective thermal conductivity constants based on mixture rules
% k1: longitudinal thermal conductivity (linear rule of mixtures)
% k2: transverse thermal conductivity (inverse rule of mixtures)
k1 = kf1*f + km1*(1-f);
k2 = 1./(f/kf2 + (1-f)/km2);

% Assemble the local thermal conductivity matrix K0.
K0(1,1,:) = k1;
K0(2,2,:) = k2;
K0(3,3,:) = k2;

% Derivatives of the effective thermal conductivity with respect to f
dK0df(1,1,:) = kf1 - km1;
dK0df(2,2,:) = -k2.^2.*(1/kf2 - 1/km2);
dK0df(3,3,:) = -k2.^2.*(1/kf2 - 1/km2);

% =========================================================================
% THERMAL EXPANSION PROPERTIES
% =========================================================================

% Effective thermal expansion constants
% a1: longitudinal thermal expansion coefficient (Schapery's rule)
% a2: transverse thermal expansion coefficient
a1 = (af1*Ef1*f + am*Em*(1-f)) ./ (Ef1*f + Em*(1-f));
a2 = (1+vf23)*af2*f + (1+vm)*am*(1-f) - a1.*v12;

% Assemble the local thermal expansion matrix A0.
A0(1,1,:) = a1;
A0(2,2,:) = a2;
A0(3,3,:) = a2;

% Derivatives of the effective thermal expansion with respect to f
da1df = ((Ef1*af1 - Em*am).*(Ef1*f + Em*(1-f)) - (af1*Ef1*f + am*Em*(1-f)).*(Ef1 - Em)) ./ (Ef1*f + Em*(1-f)).^2;
da2df = (1+vf23)*af2 - (1+vm)*am - (da1df.*v12 + a1.*dv12);

% Assemble derivative of the thermal expansion matrix with respect to f
dA0df(1,1,:) = da1df;
dA0df(2,2,:) = da2df;
dA0df(3,3,:) = da2df;

end