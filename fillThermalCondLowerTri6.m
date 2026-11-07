function Ka = fillThermalCondLowerTri6(C22, g, l, m, n)
% Construct symmetric 3x3 tensor as a 6xN matrix (Voigt order)
l = l(:).'; m = m(:).'; n = n(:).';
C22 = C22(:).'; g = g(:).';
nElem = numel(l);
Ka = zeros(6, nElem);
Ka(1,:) = C22 + g.*l.^2;   % (1,1)
Ka(2,:) = g.*l.*m;        % (2,1)
Ka(3,:) = g.*l.*n;        % (3,1)
Ka(4,:) = C22 + g.*m.^2;   % (2,2)
Ka(5,:) = g.*m.*n;        % (3,2)
Ka(6,:) = C22 + g.*n.^2;   % (3,3)
end