function Kth = H8ThermalKth(lx, ly, lz, kappa)
% H8ThermalKth: Compute the thermal conductivity matrix of a 3D eight-node hexahedral element.
% Inputs:
% lx, ly, lz: Element dimensions in the x, y, and z directions.
% kappa: 3x3 thermal conductivity tensor in the global coordinate system.

% Element half-lengths.
a = lx/2; b = ly/2; c = lz/2;

% Gauss points and weights for 2x2x2 integration.
xx=[-1/sqrt(3), 1/sqrt(3)]; yy = xx;zz=xx;
ww=[1,1];

% Initialize the element matrix.
Kth = zeros(8,8);   % Element thermal conductivity matrix.
detJ=a*b*c;
% Natural coordinates of the eight nodes.
xic   = [-1  1  1 -1 -1  1  1 -1];
etac  = [-1 -1  1  1 -1 -1  1  1];
zetac = [-1 -1 -1 -1  1  1  1  1];

%% Compute Kth by volume integration
for ii=1:length(xx)
    for jj=1:length(yy)
        for kk=1:length(zz)
            xi = xx(ii);eta = yy(jj);zeta = zz(kk);
            dNx=xic.*(1+eta*etac).*(1+zeta*zetac)/8/a;
            dNy=(1+xi*xic).*etac.*(1+zeta*zetac)/8/b;
            dNz=(1+xi*xic).*(1+eta*etac).*zetac/8/c;
            
            % Assemble the thermal gradient matrix Bth (3x8).
            % Temperature gradient: g = Bth * Te.
            Bth = [dNx; dNy; dNz]; 
            weight = ww(ii)*ww(jj)*ww(kk)*detJ;
            
            % Accumulate the element conductivity matrix.
            Kth = Kth + weight * (Bth' * kappa * Bth);
        end
    end
end
end