function [xp,dxp]=MaterialInterpolationKappa(x)
eps = 1e-3;
xp = eps + (1-eps)*x;
dxp = 1-eps;
end