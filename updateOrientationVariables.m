%% -------------------------------------- DecayedPMove local-asymptote orientation update
function znew = updateOrientationVariables(dfdz,z,Iter)
% Updates one orientation component vector z = px, py or pz.
% PMove type: DecayedPMove, move = (1/8)*0.95^(Iter-1).
% Asymptote type: LocalAsymptote, i.e. L and U are reset around current z.

zMin = -1.0;
zMax =  1.0;
move = 1/8*0.95^(Iter-1);
AsymInit = 0.5;
albefa = 0.1;

% Compute local/reset asymptotes L and U.
L = z - AsymInit*(zMax-zMin);
U = z + AsymInit*(zMax-zMin);

% Compute move-limited bounds alpha and beta.
zzz1 = L + albefa*(z-L);
zzz2 = z - move*(zMax-zMin);
alpha = max(zMin,max(zzz1,zzz2));

zzz1 = U - albefa*(U-z);
zzz2 = z + move*(zMax-zMin);
beta = min(zMax,min(zzz1,zzz2));

% Solve the separable approximate subproblem.
p = (U-z).^2.*max(dfdz,0);
q = (z-L).^2.*(-min(dfdz,0));
p = sqrt(p);
q = sqrt(q);

% Numerical protection for zero-gradient entries where p + q = 0.
den = p + q;
zCnd = z;
idx = den > 1e-12;
zCnd(idx) = (L(idx).*p(idx) + U(idx).*q(idx))./den(idx);

znew = max(alpha,min(beta,zCnd));
end