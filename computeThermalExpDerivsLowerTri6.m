function [dAl, dAm, dAn] = computeThermalExpDerivsLowerTri6(h, l, m, n)
% Directional derivatives of thermal expansion vector (6xN components)
l = l(:).'; m = m(:).'; n = n(:).'; h = h(:).';
nElem = numel(l);
dAl = zeros(6, nElem);
dAm = zeros(6, nElem);
dAn = zeros(6, nElem);

% Derivatives w.r.t l
dAl(1,:) = 2*h.*l;
dAl(5,:) = 2*h.*n;
dAl(6,:) = 2*h.*m;

% Derivatives w.r.t m
dAm(2,:) = 2*h.*m;
dAm(4,:) = 2*h.*n;
dAm(6,:) = 2*h.*l;

% Derivatives w.r.t n
dAn(3,:) = 2*h.*n;
dAn(4,:) = 2*h.*m;
dAn(5,:) = 2*h.*l;
end