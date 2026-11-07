function [Pe] = H8ThermalLoadVector(lx, ly, lz, Q_source, q_vec)
% H8ThermalLoadVector: Compute the equivalent thermal load vector of a 3D eight-node hexahedral element.
% Inputs:
% lx, ly, lz : Element dimensions in the x, y, and z directions.
% Q_source   : Volumetric heat-generation rate (scalar, W/m^3).
% q_vec      : Six-component surface heat-flux vector [q1, q2, q3, q4, q5, q6] (W/m^2).
%              Components correspond to the 1:x-, 2:x+, 3:y-, 4:y+, 5:z-, and 6:z+ faces.

Pe = zeros(8,1);
a = lx/2; b = ly/2; c = lz/2;
detJ = a*b*c;

% Gauss points and weights for two-point integration.
gp = [-1/sqrt(3), 1/sqrt(3)];
gw = [1, 1];

% Natural coordinates of the eight nodes.
xic   = [-1  1  1 -1 -1  1  1 -1];
etac  = [-1 -1  1  1 -1 -1  1  1];
zetac = [-1 -1 -1 -1  1  1  1  1];

%% 1. Volumetric heat-source contribution
if Q_source ~= 0
    for ii=1:2
        for jj=1:2
            for kk=1:2
                xi = gp(ii); eta = gp(jj); zeta = gp(kk);
                Nth = (1+xi*xic).*(1+eta*etac).*(1+zeta*zetac)/8;
                weight = gw(ii)*gw(jj)*gw(kk) * detJ;
                Pe = Pe + weight * (Nth' * Q_source);
            end
        end
    end
end

%% 2. Surface heat-flux contribution
% Loop over the six element faces represented in q_vec.
for face_id = 1:6
    q_val = q_vec(face_id);
    
    if q_val ~= 0
        % Perform 2x2 Gauss integration over the current face.
        for i=1:2
            for j=1:2
                % Set the natural coordinates and surface Jacobian dA for the current face.
                switch face_id
                    case 1; xi=-1; eta=gp(i); zeta=gp(j); dA=b*c; % x- face
                    case 2; xi= 1; eta=gp(i); zeta=gp(j); dA=b*c; % x+ face
                    case 3; eta=-1; xi=gp(i); zeta=gp(j); dA=a*c; % y- face
                    case 4; eta= 1; xi=gp(i); zeta=gp(j); dA=a*c; % y+ face
                    case 5; zeta=-1; xi=gp(i); eta=gp(j); dA=a*b; % z- face
                    case 6; zeta= 1; xi=gp(i); eta=gp(j); dA=a*b; % z+ face
                end
                
                % Evaluate the shape functions at the surface Gauss point.
                Nth = (1+xi*xic).*(1+eta*etac).*(1+zeta*zetac)/8;
                
                % Accumulate the surface heat-flux load.
                weight_area = gw(i)*gw(j) * dA;
                Pe = Pe + weight_area * (Nth' * q_val);
            end
        end
    end
end

end
%% Example usage
% Initialize the global thermal load vector.
% P_global = zeros(total_nodes, 1);
% 
% % Loop over all elements.
% for e = 1:num_elements
%     % Define the volumetric heat source and face heat fluxes for the current element.
%     Q = ...; q = ...;
% 
%     % Compute the element thermal load vector.
%     Pe = H8ThermalLoadVector(lx, ly, lz, Q, q);
% 
%     % Assemble the element contribution into the global vector.
%     nodes = element_nodes(e, :);
%     P_global(nodes) = P_global(nodes) + Pe;