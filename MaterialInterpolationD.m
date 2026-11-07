function [xp,dxp]=MaterialInterpolationD(x)
eps = 1e-6;
q = 8;
xp = eps+(1-eps)*x./(1+q*(1-x));
dxp = (1-eps)*(1+q)/(1+q*(1-x)).^2;
end