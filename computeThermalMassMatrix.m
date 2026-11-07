%% -------------------------------------------------------------------------
function M_T = computeThermalMassMatrix(lx, ly, lz)
% Compute thermal mass matrix M_T = integral N^T * N dOmega for H8 element
% Uses 2x2x2 Gauss quadrature (same as FE analysis)
a = lx/2; b = ly/2; c = lz/2;
gp = [-1/sqrt(3), 1/sqrt(3)];
gw = [1, 1];
xic   = [-1  1  1 -1 -1  1  1 -1];
etac  = [-1 -1  1  1 -1 -1  1  1];
zetac = [-1 -1 -1 -1  1  1  1  1];
M_T = zeros(8,8);
for ii = 1:2
    for jj = 1:2
        for kk = 1:2
            xi = gp(ii); eta = gp(jj); zeta = gp(kk);
            N = (1+xi*xic).*(1+eta*etac).*(1+zeta*zetac)/8;
            weight = gw(ii)*gw(jj)*gw(kk) * (a*b*c);
            M_T = M_T + weight * (N' * N);
        end
    end
end
end