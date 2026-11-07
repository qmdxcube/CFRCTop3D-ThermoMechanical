%% FUNCTION mgcg - MULTIGRID PRECONDITIONED CONJUGATE GRADIENTS
function [i,relres,u] = mgcg(A,b,u,Lfac,Ufac,Pu,nl,nswp,tol,maxiter)
r = b - A{1,1}*u;
normb = norm(b);
normr = norm(r);

% Robust handling of a zero initial residual / zero right-hand side.
if normr == 0
    i = 0;
    relres = 0;
    return;
end

if normb == 0
    res0 = normr;
else
    res0 = normb;
end

% Jacobi smoother
omega = 0.2;
invD = cell(nl-1,1);
for l = 1:nl-1
    invD{l,1}= 1./spdiags(A{l,1},0);
end
for i = 1:1e6
    z = VCycle(A,r,Lfac,Ufac,Pu,1,nl,invD,omega,nswp);
    rho = r'*z;
    if i == 1
        p = z;
    else
        beta = rho/rho_p;
        p = beta*p + z;
    end
    q = A{1,1}*p;
    dpr = p'*q;
    alpha = rho/dpr;
    u = u+alpha*p;
    r = r-alpha*q;
    rho_p = rho;
    relres = norm(r)/res0;
    if relres < tol || i>=maxiter
        break
    end
end
end