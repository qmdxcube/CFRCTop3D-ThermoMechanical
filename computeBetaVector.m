function Beta = computeBetaVector(C, A)
% computeBetaVector
% -------------------------------------------------------------------------
% Efficiently computes the stress-temperature coefficient vector (Beta = D * A)
% using the 21 lower-triangular stiffness components (C) and the 6 components
% of the thermal expansion vector (A) for all elements.
%
% Inputs:
%   C - 21 lower-triangular components of stiffness matrix, size: 21 x nElem
%   A - 6 components of thermal expansion vector, size: 6 x nElem
%
% Output:
%   Beta - Calculated 6-component stress-temperature vector, size: 6 x nElem
% -------------------------------------------------------------------------

nElem = size(C, 2);
Beta = zeros(6, nElem);

% Direct mapping of matrix-vector multiplication for symmetric 6x6 matrix 
% represented as a 21-component lower triangle.
Beta(1,:) = C(1,:).*A(1,:) + C(2,:).*A(2,:)  + C(3,:).*A(3,:)  + C(4,:).*A(4,:)  + C(5,:).*A(5,:)  + C(6,:).*A(6,:);
Beta(2,:) = C(2,:).*A(1,:) + C(7,:).*A(2,:)  + C(8,:).*A(3,:)  + C(9,:).*A(4,:)  + C(10,:).*A(5,:) + C(11,:).*A(6,:);
Beta(3,:) = C(3,:).*A(1,:) + C(8,:).*A(2,:)  + C(12,:).*A(3,:) + C(13,:).*A(4,:) + C(14,:).*A(5,:) + C(15,:).*A(6,:);
Beta(4,:) = C(4,:).*A(1,:) + C(9,:).*A(2,:)  + C(13,:).*A(3,:) + C(16,:).*A(4,:) + C(17,:).*A(5,:) + C(18,:).*A(6,:);
Beta(5,:) = C(5,:).*A(1,:) + C(10,:).*A(2,:) + C(14,:).*A(3,:) + C(17,:).*A(4,:) + C(19,:).*A(5,:) + C(20,:).*A(6,:);
Beta(6,:) = C(6,:).*A(1,:) + C(11,:).*A(2,:) + C(15,:).*A(3,:) + C(18,:).*A(4,:) + C(20,:).*A(5,:) + C(21,:).*A(6,:);

end