%% FUNCTION prepcoarseKth - PREPARE MG PROLONGATION OPERATOR FOR 1-DOF (THERMAL)
function [Pu] = prepcoarseKth(nex,ney,nez)
% Assemble state variable prolongation for 1-DOF
nexc = nex/2; neyc = ney/2; nezc = nez/2;
% Each coarse-grid node can interpolate to at most 3x3x3 = 27 neighboring fine-grid nodes.
maxnum = (nexc+1)*(neyc+1)*(nezc+1)*27; 
iP = zeros(maxnum,1); jP = zeros(maxnum,1); sP = zeros(maxnum,1);

% Weights for fixed distances to neighbors on a structured grid
vals = [1,0.5,0.25,0.125];
cc = 0;
for nx = 1:nexc+1
    for ny = 1:neyc+1
        for nz = 1:nezc+1
            % Coordinate on coarse grid
            col = (nz-1)*(nexc+1)*(neyc+1)+(ny-1)*(nexc+1)+nx;
            % Coordinate on fine grid
            nx1 = nx*2 - 1; ny1 = ny*2 - 1; nz1 = nz*2 - 1;
            % Loop over fine nodes within the rectangular domain
            for k = max(nx1-1,1):min(nx1+1,nex+1)
                for l = max(ny1-1,1):min(ny1+1,ney+1)
                    for h = max(nz1-1,1):min(nz1+1,nez+1)
                        row = (h-1)*(nex+1)*(ney+1)+(l-1)*(nex+1)+k;
                        % Based on squared dist assign weights: 1.0 0.5 0.25 0.125
                        ind = 1+((nx1-k)^2+(ny1-l)^2+(nz1-h)^2);
                        
                        % Thermal problem has one DOF per node; no factor of 3 is applied to row or column indices.
                        cc = cc+1; 
                        iP(cc) = row; 
                        jP(cc) = col; 
                        sP(cc) = vals(ind);
                    end
                end
            end
        end
    end
end
% Assemble matrices
Pu = sparse(iP(1:cc),jP(1:cc),sP(1:cc));
end