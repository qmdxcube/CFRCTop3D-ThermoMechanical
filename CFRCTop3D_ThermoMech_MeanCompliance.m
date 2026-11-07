function CFRCTop3D_ThermoMech_MeanCompliance(nelx,nely,nelz,he,volfrac,contentfrac,rxmin,rymin,rpmin,cgtol,cgmax,maxiter,q,FmSum)
% Decoupled update method: [x;y] are updated first, then [px;py;pz] are updated separately.
% Thermo-mechanical coupled topology optimization of continuous fiber reinforced composites.
%% MGPCG
% Use only coarsening levels that are valid in all three mesh directions.
meshDims = [nelx, nely, nelz];
if any(meshDims <= 0) || any(meshDims ~= round(meshDims))
    error('nelx, nely, and nelz must be positive integers.');
end
numFactorsOfTwo = arrayfun(@(n) sum(factor(n) == 2), meshDims);
if min(numFactorsOfTwo) < 1
    error('nelx, nely, and nelz must all be divisible by 2 for the geometric multigrid hierarchy.');
end
nl = min(numFactorsOfTwo) + 1;

%% Material properties (Thermo-mechanical)
Em = 4.35e9; vm = 0.37; km1 = 0.19; km2 = 0.19; am = 43.92e-6; % Material properties of matrix
Ef1 = 233.13e9; Ef2 = 23.11e9; vf12 = 0.2; vf23 = 0.4; Gf12 = 8.97e9; kf1 = 8.8; kf2 = 2.0; af1 = -0.54e-6; af2 = 10.08e-6; % Material properties of fiber
ymax = 0.5;  % Upper bound for fiber volume fraction

%% Define support and mechanical load
nElem = nelx*nely*nelz; % Number of elements
nNode = (nelx+1)*(nely+1)*(nelz+1); % number of nodes
nDof = 3*nNode; % Number of dofs

% USER-DEFINED SUPPORT FIXED DOFs
[kfs,jfs,ifs] = ndgrid(1:nelz+1,1:nely+1,1); % Coordinates for fixed face (x=0)
fixednid = sub2ind([nelz+1,nely+1,nelx+1],kfs(:),jfs(:),ifs(:)); % Fixed Node IDs
fixeddofs = [3*fixednid(:)-2;3*fixednid(:)-1;3*fixednid(:)]; % Fixed DOFs

% USER-DEFINED LOAD DOFs
[kl,jl,il] = ndgrid(nelz+1,1:nely+1,1:nelx+1); % Coordinates for loaded edge
loadnid = sub2ind([nelz+1,nely+1,nelx+1],kl(:),jl(:),il(:)); % Loading Node IDs
loaddofs = 3*loadnid(:); % Loading DOFs
if nargin < 14 || isempty(FmSum)
    FmSum = -1000000;
end
if ~isscalar(FmSum)
    error('FmSum must be a scalar total mechanical load.');
end
Fm = sparse(loaddofs,1,FmSum/numel(loaddofs),nDof,1);
U = zeros(nDof,1);
accom = zeros(nNode,1);
PTOri = zeros(nNode,1);

%% Define heat flux and thermal boundary
% 6-component vector [q1, q2, q3, q4, q5, q6] (W/m^2), corresponding to: 
% 1:x-, 2:x+, 3:y-, 4:y+, 5:z-, 6:z+ face
if nargin < 13 || isempty(q)
    q = 2.0;
end
if ~isscalar(q)
    error('q must be a scalar heat flux applied on the z+ face.');
end
q1 = [0,0,0,0,0,q]; 
[kq1,jq1,iq1] = ndgrid(nelz,1:nely,1:nelx); % Elements subject to heat flux
qelem1 = sub2ind([nelz,nely,nelx],kq1(:),jq1(:),iq1(:));

BoundT = 20; % Prescribed temperature on thermal boundary nodes
Tini = 20;   % Stress-free reference temperature for thermal strain
[kBN,jBN,iBN] = ndgrid(1,1:nely+1,1:nelx+1);
BoundNode = sub2ind([nelz+1,nely+1,nelx+1],kBN(:),jBN(:),iBN(:));

%% Passive elements definition
Passive = zeros(nelz,nely,nelx);
passiveThickness = max(1, round(nelz/20));
Passive(nelz-passiveThickness+1:nelz,1:nely,1:nelx) = 2; % 1 for void, 2 for solid (non-designable domain)

% Extract active design domain indices
actv = find(Passive == 0); % Find linear indices of all designable elements
nActv = length(actv);      % Number of active variables

passiveSolidFrac = nnz(Passive == 2)/nElem;
if volfrac <= passiveSolidFrac || volfrac > 1
    error('volfrac must be in (%.4f, 1] because passive solid elements are included in the volume constraint.', passiveSolidFrac);
end

%% Finite element analysis preparation
nodenrs = int32(reshape(1:nNode,1+nelz,1+nely,1+nelx));
edofVec = reshape(nodenrs(1:end-1,1:end-1,1:end-1),nElem,1);
nids = edofVec+int32([0,(nely+1)*(nelz+1),(nely+2)*(nelz+1),nelz+1,1+[0,(nely+1)*(nelz+1),(nely+2)*(nelz+1),nelz+1]]);
edofMat = 3*repelem(nids,1,3)+repmat(int32([-2,-1,0]),1,8);

[sI, sII] = find(tril(ones(24)));
[iK,jK] = deal(edofMat(:,sI)',edofMat(:,sII)');
[sI, sII] = find(tril(ones(8)));
[iKth,jKth] = deal(nids(:,sI)',nids(:,sII)');
iFth = double(reshape(edofMat',[],1));
[sI, sII] = find((ones(24,8)));
[iH,jH] = deal(edofMat(:,sI)',nids(:,sII)');

% Preparation for stiffness matrix templates
numD = 21; % Number of independent components
KET = zeros(24,24,numD); % Template Stiffness Matrices (TSMs) in matrix form
KETc = zeros(300, numD); % Template Stiffness Matrices (TSMs) in vector form
trilIndex = tril(true(24));
[rowD, colD] = find(tril(ones(6)));
for n = 1:numD
    DT = zeros(6,6);
    i = rowD(n); j = colD(n);
    DT(i,j) = 1; DT(j,i) = 1;
    KEt = H8BasicKe(he,he,he,DT);
    KET(:,:,n) = KEt;
    KETc(:,n) = KEt(trilIndex);
end

% Preparation for thermal conductivity matrix templates
numKa = 6; % Number of independent components
KthET = zeros(8,8,numKa);
KthETc = zeros(36, numKa); 
trilIndexKth = tril(true(8));
[rowKa, colKa] = find(tril(ones(3)));
for n = 1:numKa
    KaT = zeros(3,3);
    i = rowKa(n); j = colKa(n);
    KaT(i,j) = 1; KaT(j,i) = 1;
    KthEt = H8ThermalKth(he,he,he,KaT);
    KthET(:,:,n) = KthEt;
    KthETc(:,n) = KthEt(trilIndexKth);
end

% Preparation for thermal-stress coupling vector templates
numH = 6; % Number of independent components
HET = zeros(24,8,numH);
HETc = zeros(24*8, numH); 
for n = 1:numH
    betaT = zeros(6,1);
    betaT(n) = 1;
    HEt = H8ThermalHth(he,he,he,betaT);
    HET(:,:,n) = HEt;
    HETc(:,n) = HEt(:);
end
HET2D = reshape(HET,24,8*numH);

% Prolongation operators (mechanical)
PuK = cell(nl-1,1);
for l = 1:nl-1
    [PuK{l,1}] = prepcoarse(nelz/2^(l-1),nely/2^(l-1),nelx/2^(l-1));
end    
NK = ones(nDof,1); NK(fixeddofs) = 0; NullK = spdiags(NK,0,nDof,nDof);

% Prolongation operators (thermal)
PuKth = cell(nl-1,1);
for l = 1:nl-1
    [PuKth{l,1}] = prepcoarseKth(nelz/2^(l-1),nely/2^(l-1),nelx/2^(l-1));
end    
NKth = ones(nNode,1); NKth(BoundNode) = 0; NullKth = spdiags(NKth,0,nNode,nNode);

%% Compute thermal load vector P array
Pe0 = H8ThermalLoadVector(he, he, he, 0, q1);
for e = qelem1'
    nodes = nids(e, :);
    PTOri(nodes) = PTOri(nodes) + Pe0;
end

%% Preparation for filtering and projection
[Hx,Hxs] = PreFILTER3D(nelx,nely,nelz,rxmin);
[Hy,Hys] = PreFILTER3D(nelx,nely,nelz,rymin);
[Hp,Hps] = PreFILTER3D(nelx,nely,nelz,rpmin);
filt = @(v,H,Hs)convn(v,H,'same')./Hs; % filtering for design variables
ifiltV = @(v,H,Hs) reshape(convn(reshape(v,nelz,nely,nelx)./Hs,H,'same'),[],1); % adjoint operation for sensitivities
prj = @(v,eta,beta) (tanh(beta*eta)+tanh(beta*(v-eta)))./(tanh(beta*eta)+tanh(beta*(1-eta))); % projection
dprj = @(v,eta,beta) beta*(1-tanh(beta*(v-eta)).^2)./(tanh(beta*eta)+tanh(beta*(1-eta)));% proj. x-derivative

%% Initialization of topology, fiber volume fraction and fiber orientation
yinit = contentfrac/volfrac;
if yinit < 0 || yinit > ymax
    error('Initial fiber volume fraction contentfrac/volfrac must lie in [0, %.3f].', ymax);
end
xinitActv = (volfrac*nElem - nnz(Passive == 2))/nActv;
x = xinitActv*ones(nelz,nely,nelx); x(Passive==1)=0; x(Passive==2)=1;
y = yinit*ones(nelz,nely,nelx);
px = ones(nelz,nely,nelx);
py = zeros(nelz,nely,nelx);
pz = zeros(nelz,nely,nelx);

%% Parameters
% Optimization parameters
iter = 0;
change = 1;
tol = 0.01;

% Heaviside projection parameters
loopbeta = 0;
heavbeta = 1;
heaveta = 0.5;

% MMA updates density on active elements and fiber fraction on all elements.
Nc = 2;
Nd = nActv + nElem; 
aa0 = 1;
aa = zeros(Nc,1);
cc = 1e4*ones(Nc,1);
dd = ones(Nc,1);

xyval = [x(actv); y(:)]; 
xyold1 = xyval;
xyold2 = xyval;
xymin = [zeros(nActv,1); zeros(nElem,1)];
xymax = [ones(nActv,1); ymax*ones(nElem,1)];
xyL = xymin;
xyU = xymax;

compliance = zeros(maxiter,1);
volumefrac = zeros(maxiter,1);
contentfracHis = zeros(maxiter,1);
OBJHis = zeros(maxiter,1);
timeHis = zeros(maxiter,8);

% MGCG history for every optimization iteration.
% Columns:
%   1 - thermal state equation
%   2 - mechanical equilibrium equation
%   3 - thermal adjoint equation
mgcgIterHis = zeros(maxiter,3);
mgcgTimeHis = zeros(maxiter,3);
mgcgCallLabels = { ...
    'Thermal state', ...
    'Mechanical equilibrium', ...
    'Thermal adjoint'};

TAug = zeros(nNode,1);
TAug(BoundNode, 1) = BoundT;
T = TAug;

%% START ITERATION
while (change>tol && iter<maxiter)
    iter = iter+1;
    loopbeta = loopbeta+1;
    
    %% FILTERING AND PROJECTION/NORMALIZATION
    tic
    xTilde = filt(x,Hx,Hxs); yPhys = filt(y,Hy,Hys);
    xPhys = prj(xTilde,heaveta,heavbeta);
    xPhys(Passive==1) = 0; % Enforce passive void elements
    xPhys(Passive==2) = 1; % Enforce passive solid elements
    
    pxTilde = filt(px,Hp,Hps); pyTilde = filt(py,Hp,Hps); pzTilde = filt(pz,Hp,Hps);
    Lp = sqrt(pxTilde.^2 + pyTilde.^2 + pzTilde.^2);
    Lp = max(Lp, 1e-12);
    pxPhys = pxTilde./Lp; pyPhys = pyTilde./Lp; pzPhys = pzTilde./Lp;
    
    timeHis(iter,1) = toc;
    
    %% FE-ANALYSIS
    tic
    l1 = pxPhys(:); m1 = pyPhys(:); n1 = pzPhys(:);
    [xmpD,dxmpD] = MaterialInterpolationD(xPhys); 
    [xmpKa,dxmpKa] = MaterialInterpolationKappa(xPhys); 
    [xmpBeta,dxmpBeta] = MaterialInterpolationBeta(xPhys); 
    
    [D0,dD0y,~,Ka0,A0,dKa0dy,dA0dy] = compositeMixingParams(km1,km2,am,Em,vm,kf1,kf2,af1,af2,Ef1,Ef2,vf12,vf23,Gf12,yPhys);
    [Dcomp, Kacomp, Acomp, dDlcomp, dDmcomp, dDncomp, dDycomp, ...
     dKaycomp, dKalcomp, dKamcomp, dKancomp, dAycomp, dAlcomp, dAmcomp, dAncomp] = ...
     rotateMaterialPropsThermoMech(D0, Ka0, A0, l1, m1, n1, dD0y, dKa0dy, dA0dy);
    
    Beta = computeBetaVector(Dcomp, Acomp);
    dBetay = computeBetaVector(dDycomp, Acomp) + computeBetaVector(Dcomp, dAycomp);
    dBetal = computeBetaVector(dDlcomp, Acomp) + computeBetaVector(Dcomp, dAlcomp);
    dBetam = computeBetaVector(dDmcomp, Acomp) + computeBetaVector(Dcomp, dAmcomp);
    dBetan = computeBetaVector(dDncomp, Acomp) + computeBetaVector(Dcomp, dAncomp);
    
    timeHis(iter,2) = toc;
    
    %% ASSEMBLE GLOBAL MATRICES
    tic
    sK = KETc * (Dcomp .* xmpD(:).');
    sKth = KthETc * (Kacomp .* xmpKa(:).');
    timeHis(iter,3) = toc;
    
    %% SOLVE THERMAL SYSTEM
    thermalStageTimer = tic;
    Kth = cell(nl,1);
    Kth{1,1} = sparse(iKth,jKth,sKth(:),nNode,nNode);
    Kth{1,1} = Kth{1,1}+Kth{1,1}'-diag(diag(Kth{1,1}));
    
    PT = PTOri;
    PT = PT - Kth{1,1}*TAug;
    PT(BoundNode) = BoundT;
    
    Kth{1,1} = NullKth'*Kth{1,1}*NullKth - (NullKth-speye(nNode,nNode));
    for l = 1:nl-1
        Kth{l+1,1} = PuKth{l,1}'*(Kth{l,1}*PuKth{l,1});
    end
    Lfacth = chol(Kth{nl,1},'lower'); Ufacth = Lfacth';
    
    if q ~= 0 || BoundT ~= 0
        mgcgTimer = tic;
        [mgcgIterHis(iter,1), ~, T] = ...
            mgcg(Kth, PT, T, Lfacth, Ufacth, PuKth, nl, 1, cgtol, cgmax);
        mgcgTimeHis(iter,1) = toc(mgcgTimer);
    else
        T = zeros(nNode,1);
        mgcgIterHis(iter,1) = 0;
        mgcgTimeHis(iter,1) = 0;
    end
    timeHis(iter,4) = toc(thermalStageTimer);
    
    %% SOLVE MECHANICAL SYSTEM WITH THERMAL LOAD
    mechanicalStageTimer = tic;
    deltaT = T - Tini;
    deltaTe = deltaT(nids)'; % 8 x nElem matrix
    
    % Assemble thermal load vector Fth using template thermo-mechanical matrices
    scaledBeta = Beta .* xmpBeta(:).';
    V3D = reshape(deltaTe, 8, 1, []) .* reshape(scaledBeta, 1, 6, []);
    sFth = HET2D * reshape(V3D, 48, []);
    sH = HETc * (Beta .* xmpBeta(:).');
    Fth = sparse(iFth,1,sFth(:),nDof,1);
    
    F = Fm + Fth;
    F = NullK * F;
    
    K = cell(nl,1);
    K{1,1} = sparse(iK,jK,sK(:),nDof,nDof);
    K{1,1} = K{1,1}+K{1,1}'-diag(diag(K{1,1}));
    K{1,1} = NullK'*K{1,1}*NullK - (NullK-speye(nDof,nDof));
    
    for l = 1:nl-1
        K{l+1,1} = PuK{l,1}'*(K{l,1}*PuK{l,1});
    end
    Lfac = chol(K{nl,1},'lower'); Ufac = Lfac';
    
    mgcgTimer = tic;
    [mgcgIterHis(iter,2), ~, U] = ...
        mgcg(K, F, U, Lfac, Ufac, PuK, nl, 1, cgtol, cgmax);
    mgcgTimeHis(iter,2) = toc(mgcgTimer);
    
    timeHis(iter,5) = toc(mechanicalStageTimer);
    
    compliance(iter,1) = F'*U;
    OBJHis(iter,1) = compliance(iter,1);
    volumefrac(iter,1) = mean(xPhys(:));
    contentfracHis(iter,1) = mean(xPhys(:).*yPhys(:));
    
    %% SENSITIVITY ANALYSIS
    sensitivityStageTimer = tic;
    % Adjoint variable solver for thermal system
    Hth = sparse(double(iH),double(jH),sH(:));
    b = Hth'*U;
    b(BoundNode) = 0;
    mgcgTimer = tic;
    [mgcgIterHis(iter,3), ~, accom] = ...
        mgcg(Kth, b, accom, Lfacth, Ufacth, PuKth, nl, 1, cgtol, cgmax);
    mgcgTimeHis(iter,3) = toc(mgcgTimer);
    
    if mod(iter,20) == 1
        scale = 1/compliance(iter,1); 
    end
    
    Ud = U(edofMat);
    
    % Mechanical sensitivity contributions
    uku = zeros(nElem, 21);
    for n = 1:21
        uku(:,n) = -sum((Ud * KET(:,:,n)) .* Ud, 2);
    end

    cex = dxmpD(:) .* sum(Dcomp.' .* uku, 2);
    cey = xmpD(:)  .* sum(dDycomp.' .* uku, 2);
    cepx = xmpD(:) .* sum(dDlcomp.' .* uku, 2);
    cepy = xmpD(:) .* sum(dDmcomp.' .* uku, 2);
    cepz = xmpD(:) .* sum(dDncomp.' .* uku, 2);
    
    % Thermal sensitivity contributions
    accomN = accom(nids);
    TN = T(nids);
    deltaTN = deltaT(nids);
    
    akthT = zeros(nElem, 6);
    uFth = zeros(nElem, 6);
    for n = 1:6
        akthT(:,n) = -2*sum((accomN * KthET(:,:,n)) .* TN, 2);
        uFth(:,n) = 2*sum((Ud * HET(:,:,n)) .* deltaTN, 2);
    end
    
    cex = cex + dxmpKa(:) .* sum(Kacomp.' .* akthT, 2) + dxmpBeta(:) .* sum(Beta.' .* uFth, 2);
    cey = cey + xmpKa(:)  .* sum(dKaycomp.' .* akthT, 2) + xmpBeta(:)  .* sum(dBetay.' .* uFth, 2);
    cepx = cepx + xmpKa(:) .* sum(dKalcomp.' .* akthT, 2) + xmpBeta(:) .* sum(dBetal.' .* uFth, 2);
    cepy = cepy + xmpKa(:) .* sum(dKamcomp.' .* akthT, 2) + xmpBeta(:) .* sum(dBetam.' .* uFth, 2);
    cepz = cepz + xmpKa(:) .* sum(dKancomp.' .* akthT, 2) + xmpBeta(:) .* sum(dBetan.' .* uFth, 2);
    
    cex = scale*cex; cey = scale*cey;
    cepx = scale*cepx; cepy = scale*cepy; cepz = scale*cepz;
    
    dcdxPhysV = cex; dcdyPhysV = cey;
    dcdpxPhysV = cepx; dcdpyPhysV = cepy; dcdpzPhysV = cepz;
    dcdxPhysV(Passive(:) ~= 0) = 0;
    
    %% FILTERING OF SENSITIVITIES
    % Objective function
    dx = dprj(xTilde(:), heaveta, heavbeta);
    dcdx = ifiltV(dcdxPhysV .* dx, Hx, Hxs);
    dcdy = ifiltV(dcdyPhysV, Hy, Hys);
    
    G = [dcdpxPhysV, dcdpyPhysV, dcdpzPhysV]; 
    P = [l1,m1,n1]; 
    dcdpTiles = (G-P.*sum(P.*G,2))./Lp(:);
    dcdpx = ifiltV(dcdpTiles(:,1), Hp, Hps);
    dcdpy = ifiltV(dcdpTiles(:,2), Hp, Hps);
    dcdpz = ifiltV(dcdpTiles(:,3), Hp, Hps);
    
    f0val = scale*compliance(iter);
    df0dxy = [dcdx(actv); dcdy];
    
    % Constraint functions
    dgvdxPhys = ones(nelz,nely,nelx); dgvdxPhys(Passive~=0)=0;
    dgvdx = ifiltV(dgvdxPhys(:).*dx/nElem, Hx, Hxs);
    dgvdy = zeros(nElem,1);
    
    dgfdxPhysV = yPhys(:); dgfdyPhysV = xPhys(:);
    dgfdxPhysV(Passive(:) ~= 0) = 0;
    dgfdx = ifiltV(dgfdxPhysV.*dx/nElem, Hx, Hxs);
    dgfdy = ifiltV(dgfdyPhysV/nElem, Hy, Hys);
    
    fval = zeros(Nc,1);
    fval(1,1) = mean(xPhys(:))/volfrac - 1;
    fval(2,1) = mean(xPhys(:).*yPhys(:))/contentfrac - 1;
    
    dfvaldxy = zeros(Nc,Nd);
    dfvaldxy(1,:) = [dgvdx(actv)', dgvdy'] / volfrac;
    dfvaldxy(2,:) = [dgfdx(actv)', dgfdy'] / contentfrac;
    
    timeHis(iter,6) = toc(sensitivityStageTimer);
    
    %% DESIGN VARIABLE UPDATE
    % Update of active density variables and fiber volume fractions
    tic
    [xymma,~,~,~,~,~,~,~,~,xyL,xyU] = mmasub(Nc,Nd,iter,xyval,xymin,xymax,xyold1,xyold2,f0val,df0dxy,fval,dfvaldxy,xyL,xyU,aa0,aa,cc,dd);
    
    x(actv) = xymma(1:nActv, 1); 
    y = reshape(xymma(nActv+1:end, 1),nelz,nely,nelx);
    
    xyold2 = xyold1; xyold1 = xyval; 
    xyval = [x(actv); y(:)];
    
    change = max(abs(xyval-xyold1));
    timeHis(iter,7) = toc;
    
    % Update of global fiber orientation variables
    tic
    df0dp = [dcdpx; dcdpy; dcdpz];
    pval = [px(:); py(:); pz(:)];
    pval = updateOrientationVariables(df0dp,pval,iter);
    
    px = reshape(pval(1:nElem),nelz,nely,nelx);
    py = reshape(pval(nElem+(1:nElem)),nelz,nely,nelx);
    pz = reshape(pval(2*nElem+(1:nElem)),nelz,nely,nelx);
    timeHis(iter,8) = toc;
    
    %% PRINT RESULTS
    fprintf([' It.:%5i Obj.:%11.4f Com.:%11.4f Vol.:%7.3f ', ...
             'fVf.:%7.3f ch.:%7.3f MGCGit.:%5i MGCGt.:%8.4f s\n'], ...
        iter, OBJHis(iter,1), compliance(iter,1), volumefrac(iter,1), ...
        contentfracHis(iter,1), change, ...
        sum(mgcgIterHis(iter,:)), sum(mgcgTimeHis(iter,:)));
    
    %% UPDATE HEAVISIDE REGULARIZATION PARAMETER
    if  heavbeta < 128 && (loopbeta >= 30 || change < tol)
        heavbeta = 2*heavbeta;
        loopbeta = 0;
        change = 1;
        fprintf('Parameter heavbeta increased to %g.\n', heavbeta);
    end
end

compliance = compliance(1:iter,1);
volumefrac = volumefrac(1:iter,1);
contentfracHis = contentfracHis(1:iter,1);
timeHis = timeHis(1:iter,:);

mgcgIterHis = mgcgIterHis(1:iter,:);
mgcgTimeHis = mgcgTimeHis(1:iter,:);
mgcgTotalIterPerOptIter = sum(mgcgIterHis,2);
mgcgTotalTimePerOptIter = sum(mgcgTimeHis,2);

%% VISUALIZATION OF THE OPTIMIZED DESIGN
% plotOptimizedDesign(he,he,he,xPhys,yPhys,pxPhys,pyPhys,pzPhys,U,T)
% plotHistory(compliance,volumefrac,contentfracHis,'Mean Compliance')

%% VALIDATION OF THE OPTIMALITY OF FIBER ORIENTATION
[~,~,S0] = compositeMixingParams(km1,km2,am,Em,vm,kf1,kf2,af1,af2,Ef1,Ef2,vf12,vf23,Gf12,yPhys);
% Optionally call extra stress evaluation if implemented
% [strain, stress, FIEle] = StrainStressCalAnisotropy(he, he, he, U, D, xPhys, pxPhys, pyPhys, pzPhys);
% Reconstruct the full 6-by-6-by-nElem stiffness matrix only for validation and saving.
D = zeros(6,6,nElem);
for k = 1:numD
    Dk = reshape(Dcomp(k,:),1,1,nElem);
    D(rowD(k),colD(k),:) = Dk;
    D(colD(k),rowD(k),:) = Dk;
end
% verifyOptimality3D(he,he,he,U,D,S0,xPhys,yPhys,pxPhys,pyPhys,pzPhys);
filename = sprintf('CFRCTop3D_ThermoMech_MGPCG_MGCGHistory_%dx%dx%d_vf%.1fcf%.2f_q%.1fFm%.1f.mat',nelx,nely,nelz,volfrac,contentfrac,q,FmSum);
save(filename, ...
    'compliance','volumefrac','contentfracHis', ...
    'xPhys','yPhys','pxPhys','pyPhys','pzPhys', ...
    'U','D','S0','he','T','q','FmSum','timeHis', ...
    'mgcgIterHis','mgcgTimeHis','mgcgCallLabels', ...
    'mgcgTotalIterPerOptIter','mgcgTotalTimePerOptIter');
