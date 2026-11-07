function A = fillThermalExpLowerTri6(C22, h, l, m, n)
% Construct rotated thermal expansion vector with Voigt engineering shear terms (2*alpha_ij)
l = l(:).'; m = m(:).'; n = n(:).';
C22 = C22(:).'; h = h(:).';
nElem = numel(l);
A = zeros(6, nElem);
A(1,:) = C22 + h.*l.^2;   % (1,1)
A(2,:) = C22 + h.*m.^2;   % (2,2)
A(3,:) = C22 + h.*n.^2;   % (3,3)
A(4,:) = 2*h.*m.*n;       % 2*(2,3)
A(5,:) = 2*h.*l.*n;       % 2*(1,3)
A(6,:) = 2*h.*l.*m;       % 2*(1,2)
end