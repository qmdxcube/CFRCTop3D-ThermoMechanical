function [dKal, dKam, dKan] = computeThermalCondDerivsLowerTri6(g, l, m, n)
% Directional derivatives of thermal conductivity (6xN components)
l = l(:).'; m = m(:).'; n = n(:).'; g = g(:).';
nElem = numel(l);
dKal = zeros(6, nElem);
dKam = zeros(6, nElem);
dKan = zeros(6, nElem);

% Derivatives w.r.t l
dKal(1,:) = 2*g.*l;
dKal(2,:) = g.*m;
dKal(3,:) = g.*n;

% Derivatives w.r.t m
dKam(2,:) = g.*l;
dKam(4,:) = 2*g.*m;
dKam(5,:) = g.*n;

% Derivatives w.r.t n
dKan(3,:) = g.*l;
dKan(5,:) = g.*m;
dKan(6,:) = 2*g.*n;
end