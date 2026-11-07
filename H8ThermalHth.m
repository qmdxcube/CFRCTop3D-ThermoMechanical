function H = H8ThermalHth(lx, ly, lz, beta)
    a = lx/2; b = ly/2; c = lz/2;
    gp = [-1/sqrt(3), 1/sqrt(3)]; gw = [1, 1];
    H = zeros(24,8);
    
    xic=  [-1 1 1 -1 -1 1 1 -1];
    etac= [-1 -1 1 1 -1 -1 1 1];
    zetac=[-1 -1 -1 -1 1 1 1 1];
    B=zeros(6,24);
    for ii=1:2
        for jj=1:2
            for kk=1:2
                xi=gp(ii); eta=gp(jj); zeta=gp(kk);
                % 1x8 shape-function vector for temperature interpolation.
                Nth = (1+xi*xic).*(1+eta*etac).*(1+zeta*zetac)/8;
                % Shape-function derivatives used to assemble the 6x24 strain-displacement matrix.
                dNx = xic.*(1+eta*etac).*(1+zeta*zetac)/8/a;
                dNy = (1+xi*xic).*etac.*(1+zeta*zetac)/8/b;
                dNz = (1+xi*xic).*(1+eta*etac).*zetac/8/c;
                
                % Assemble the strain-displacement matrix B (6x24).
                B(1,1:3:24)=dNx;B(2,2:3:24)=dNy;B(3,3:3:24)=dNz;
                B(4,2:3:24)=dNz;B(4,3:3:24)=dNy;
                B(5,1:3:24)=dNz;B(5,3:3:24)=dNx;
                B(6,1:3:24)=dNy;B(6,2:3:24)=dNx;
                
                weight = gw(ii)*gw(jj)*gw(kk) * (a*b*c);
                % Accumulate the thermo-mechanical coupling template contribution.
                % H = integral B^T * beta * N dOmega, each column multiplies by shape function
                H = H + weight * ((B' * beta) * Nth);
            end
        end
    end
end