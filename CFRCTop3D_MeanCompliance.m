function CFRCTop3D_MeanCompliance(nelx,nely,nelz,he,volfrac,contentfrac,rxmin,rymin,rpmin,cgtol,cgmax,maxiter,FmSum)
% CFRCTop3D_MeanCompliance
% -------------------------------------------------------------------------
% Mechanical-load-only topology optimization of continuous-fiber-reinforced
% composites using:
%   1) density x,
%   2) fiber volume fraction y, and
%   3) fiber orientation components px, py, pz.
%
% Objective:
%   Minimum mean compliance, C_m = F''*U.
%
% The thermal conduction problem, thermal expansion, thermal load vector,
% thermal adjoint equation, and all thermal sensitivity terms are removed.
%
% Example:
% CFRCTop3D_MeanCompliance(80,40,40,1/80,0.30,0.09,1.5,1.5,1.5,1e-5,1000,200,-1e6)
%
% Required external functions:
%   PreFILTER3D
%   MaterialInterpolationD
%   computeConstitutiveMatrix3D
%   ComputeElasticityMatrices3DLowerTri21
%   H8BasicKe
%   prepcoarse
%   mgcg
%   mmasub
%   updateOrientationVariables
% -------------------------------------------------------------------------

%% MGPCG hierarchy
% Use only coarsening levels that are valid in all three mesh directions.
meshDims = [nelx,nely,nelz];
if any(meshDims <= 0) || any(meshDims ~= round(meshDims))
    error('nelx, nely, and nelz must be positive integers.');
end
numFactorsOfTwo = arrayfun(@(n) sum(factor(n) == 2), meshDims);
if min(numFactorsOfTwo) < 1
    error('nelx, nely, and nelz must all be divisible by 2 for the geometric multigrid hierarchy.');
end
nl = min(numFactorsOfTwo) + 1;

%% Mechanical material properties
Em   = 4.35e9;
vm   = 0.37;

Ef1  = 233.13e9;
Ef2  = 23.11e9;
vf12 = 0.20;
vf23 = 0.40;
Gf12 = 8.97e9;

ymax = 0.5; % Upper bound for fiber volume fraction

%% Input checks
if nargin < 13 || isempty(FmSum)
    FmSum = -1.0e6;
end

if ~isscalar(FmSum) || ~isfinite(FmSum)
    error('FmSum must be a finite scalar total mechanical load.');
end

if volfrac <= 0 || volfrac > 1
    error('volfrac must lie in (0,1].');
end

if contentfrac <= 0
    error('contentfrac must be positive.');
end

if contentfrac > ymax*volfrac
    error(['contentfrac cannot exceed ymax*volfrac = %.6g. ', ...
           'Current contentfrac = %.6g.'], ymax*volfrac, contentfrac);
end

%% Define support and mechanical load
nElem = nelx*nely*nelz;
nNode = (nelx+1)*(nely+1)*(nelz+1);
nDof  = 3*nNode;

% Fixed face: x = 0
[kfs,jfs,ifs] = ndgrid(1:nelz+1,1:nely+1,1);
fixednid = sub2ind([nelz+1,nely+1,nelx+1], ...
    kfs(:),jfs(:),ifs(:));
fixeddofs = [3*fixednid(:)-2; ...
             3*fixednid(:)-1; ...
             3*fixednid(:)];

% Mechanical load applied in the z direction on the selected edge/face
[kl,jl,il] = ndgrid(nelz+1,1:nely+1,1:nelx+1);
loadnid  = sub2ind([nelz+1,nely+1,nelx+1], ...
    kl(:),jl(:),il(:));
loaddofs = 3*loadnid(:);

Fm = sparse(loaddofs,1,FmSum/numel(loaddofs),nDof,1);
U  = zeros(nDof,1);

%% Passive elements
Passive = zeros(nelz,nely,nelx);

passiveThickness = max(1,round(nelz/20));
Passive(nelz-passiveThickness+1:nelz,1:nely,1:nelx) = 2;
% Passive = 0: active design element
% Passive = 1: prescribed void
% Passive = 2: prescribed solid

actv  = find(Passive == 0);
nActv = numel(actv);

passiveSolidFrac = nnz(Passive == 2)/nElem;

if volfrac <= passiveSolidFrac || volfrac > 1
    error(['volfrac must be in (%.6f,1] because passive solid elements ', ...
           'are included in the volume constraint.'], passiveSolidFrac);
end

%% Finite-element preparation
nodenrs = int32(reshape(1:nNode,1+nelz,1+nely,1+nelx));

edofVec = reshape( ...
    nodenrs(1:end-1,1:end-1,1:end-1), ...
    nElem,1);

nids = edofVec + int32([ ...
    0, ...
    (nely+1)*(nelz+1), ...
    (nely+2)*(nelz+1), ...
    nelz+1, ...
    1 + [ ...
        0, ...
        (nely+1)*(nelz+1), ...
        (nely+2)*(nelz+1), ...
        nelz+1]]);

edofMat = 3*repelem(nids,1,3) + ...
    repmat(int32([-2,-1,0]),1,8);

% Lower-triangular global stiffness assembly indices
[sI,sII] = find(tril(ones(24)));
[iK,jK] = deal(edofMat(:,sI)',edofMat(:,sII)');

%% Template stiffness matrices
numD = 21;

KET  = zeros(24,24,numD);
KETc = zeros(300,numD);

trilIndex = tril(true(24));
[rowD,colD] = find(tril(ones(6)));

for iTemplate = 1:numD
    DT = zeros(6,6);

    i = rowD(iTemplate);
    j = colD(iTemplate);

    DT(i,j) = 1;
    DT(j,i) = 1;

    KEt = H8BasicKe(he,he,he,DT);

    KET(:,:,iTemplate) = KEt;
    KETc(:,iTemplate)  = KEt(trilIndex);
end

%% Mechanical multigrid operators
PuK = cell(nl-1,1);

for level = 1:nl-1
    PuK{level,1} = prepcoarse( ...
        nelz/2^(level-1), ...
        nely/2^(level-1), ...
        nelx/2^(level-1));
end

NK = ones(nDof,1);
NK(fixeddofs) = 0;
NullK = spdiags(NK,0,nDof,nDof);

%% Filters and projection
[Hx,Hxs] = PreFILTER3D(nelx,nely,nelz,rxmin);
[Hy,Hys] = PreFILTER3D(nelx,nely,nelz,rymin);
[Hp,Hps] = PreFILTER3D(nelx,nely,nelz,rpmin);

filt = @(v,H,Hs) convn(v,H,'same')./Hs;

ifiltV = @(v,H,Hs) reshape( ...
    convn(reshape(v,nelz,nely,nelx)./Hs,H,'same'), ...
    [],1);

prj = @(v,eta,beta) ...
    (tanh(beta*eta)+tanh(beta*(v-eta)))./ ...
    (tanh(beta*eta)+tanh(beta*(1-eta)));

dprj = @(v,eta,beta) ...
    beta*(1-tanh(beta*(v-eta)).^2)./ ...
    (tanh(beta*eta)+tanh(beta*(1-eta)));

%% Initialize design variables
yinit = contentfrac/volfrac;

xinitActv = ...
    (volfrac*nElem - nnz(Passive == 2))/nActv;

x = xinitActv*ones(nelz,nely,nelx);
x(Passive == 1) = 0;
x(Passive == 2) = 1;

y = yinit*ones(nelz,nely,nelx);

px = ones(nelz,nely,nelx);
py = zeros(nelz,nely,nelx);
pz = zeros(nelz,nely,nelx);

%% Optimization parameters
iter    = 0;
change  = 1;
tol     = 0.01;

% Heaviside projection
loopbeta = 0;
heavbeta = 1;
heaveta  = 0.5;

%% MMA settings for x and y
Nc = 2;
Nd = nActv + nElem;

aa0 = 1;
aa  = zeros(Nc,1);
cc  = 1e4*ones(Nc,1);
dd  = ones(Nc,1);

xyval  = [x(actv);y(:)];
xyold1 = xyval;
xyold2 = xyval;

xymin = [zeros(nActv,1);zeros(nElem,1)];
xymax = [ones(nActv,1);ymax*ones(nElem,1)];

xyL = xymin;
xyU = xymax;

%% History arrays
compliance      = zeros(maxiter,1);
volumefrac      = zeros(maxiter,1);
contentfracHis  = zeros(maxiter,1);
OBJHis          = zeros(maxiter,1);

% Columns:
% 1 filtering/projection
% 2 material transformation
% 3 stiffness assembly
% 4 mechanical MGPCG solve
% 5 sensitivity analysis
% 6 MMA update of x and y
% 7 fiber-orientation update
timeHis = zeros(maxiter,7);

timeLabels = { ...
    'FilteringProjection', ...
    'MaterialTransformation', ...
    'StiffnessAssembly', ...
    'MechanicalMGPCG', ...
    'SensitivityAnalysis', ...
    'MMAUpdateXY', ...
    'OrientationUpdate'};

% Detailed MGCG history.
% The mechanical-only formulation calls MGCG once in every optimization
% iteration to solve K*U = F.
mgcgIterHis = zeros(maxiter,1);
mgcgTimeHis = zeros(maxiter,1);
mgcgCallLabels = {'Mechanical equilibrium'};

%% Optimization loop
while change > tol && iter < maxiter
    iter = iter + 1;
    loopbeta = loopbeta + 1;

    %% Filtering, projection, and orientation normalization
    tic;

    xTilde = filt(x,Hx,Hxs);
    yPhys  = filt(y,Hy,Hys);

    xPhys = prj(xTilde,heaveta,heavbeta);

    xPhys(Passive == 1) = 0;
    xPhys(Passive == 2) = 1;

    pxTilde = filt(px,Hp,Hps);
    pyTilde = filt(py,Hp,Hps);
    pzTilde = filt(pz,Hp,Hps);

    Lp = sqrt(pxTilde.^2 + pyTilde.^2 + pzTilde.^2);
    Lp = max(Lp,1e-12);

    pxPhys = pxTilde./Lp;
    pyPhys = pyTilde./Lp;
    pzPhys = pzTilde./Lp;

    timeHis(iter,1) = toc;

    %% Mechanical constitutive matrices
    tic;

    l1 = pxPhys(:);
    m1 = pyPhys(:);
    n1 = pzPhys(:);

    [xmpD,dxmpD] = MaterialInterpolationD(xPhys);

    [D0,dD0y,~] = computeConstitutiveMatrix3D( ...
        Em,vm,Ef1,Ef2,vf12,vf23,Gf12,yPhys);

    [Dcomp,dDlcomp,dDmcomp,dDncomp,dDycomp] = ...
        ComputeElasticityMatrices3DLowerTri21( ...
            D0,l1,m1,n1,dD0y);

    timeHis(iter,2) = toc;

    %% Assemble global mechanical stiffness
    tic;

    sK = KETc*(Dcomp.*xmpD(:).');

    timeHis(iter,3) = toc;

    %% Mechanical analysis under force load only
    mechanicalStageTimer = tic;

    F = NullK*Fm;

    K = cell(nl,1);

    K{1,1} = sparse(iK,jK,sK(:),nDof,nDof);
    K{1,1} = K{1,1} + K{1,1}' - diag(diag(K{1,1}));

    K{1,1} = ...
        NullK'*K{1,1}*NullK - ...
        (NullK-speye(nDof,nDof));

    for level = 1:nl-1
        K{level+1,1} = ...
            PuK{level,1}'*(K{level,1}*PuK{level,1});
    end

    Lfac = chol(K{nl,1},'lower');
    Ufac = Lfac';

    mgcgTimer = tic;
    [mgcgIterHis(iter,1),~,U] = mgcg( ...
        K,F,U,Lfac,Ufac,PuK,nl,1,cgtol,cgmax);
    mgcgTimeHis(iter,1) = toc(mgcgTimer);

    % Full mechanical-analysis time, including hierarchy construction,
    % coarse-grid factorization, and the MGCG solve.
    timeHis(iter,4) = toc(mechanicalStageTimer);

    %% Objective and constraints
    compliance(iter,1)     = F'*U;
    OBJHis(iter,1)         = compliance(iter,1);
    volumefrac(iter,1)     = mean(xPhys(:));
    contentfracHis(iter,1) = mean(xPhys(:).*yPhys(:));

    %% Mechanical sensitivity analysis
    tic;

    if mod(iter,20) == 1
        scale = 1/compliance(iter,1);
    end

    Ud = U(edofMat);

    uku = zeros(nElem,numD);

    for iTemplate = 1:numD
        uku(:,iTemplate) = ...
            -sum((Ud*KET(:,:,iTemplate)).*Ud,2);
    end

    % Purely mechanical sensitivity contributions
    cex  = dxmpD(:).*sum(Dcomp.'.*uku,2);
    cey  = xmpD(:).*sum(dDycomp.'.*uku,2);
    cepx = xmpD(:).*sum(dDlcomp.'.*uku,2);
    cepy = xmpD(:).*sum(dDmcomp.'.*uku,2);
    cepz = xmpD(:).*sum(dDncomp.'.*uku,2);

    cex  = scale*cex;
    cey  = scale*cey;
    cepx = scale*cepx;
    cepy = scale*cepy;
    cepz = scale*cepz;

    dcdxPhysV  = cex;
    dcdyPhysV  = cey;
    dcdpxPhysV = cepx;
    dcdpyPhysV = cepy;
    dcdpzPhysV = cepz;

    dcdxPhysV(Passive(:) ~= 0) = 0;

    %% Sensitivity filtering
    dx = dprj(xTilde(:),heaveta,heavbeta);

    dcdx = ifiltV(dcdxPhysV.*dx,Hx,Hxs);
    dcdy = ifiltV(dcdyPhysV,Hy,Hys);

    G = [dcdpxPhysV,dcdpyPhysV,dcdpzPhysV];
    P = [l1,m1,n1];

    dcdpTilde = ...
        (G-P.*sum(P.*G,2))./Lp(:);

    dcdpx = ifiltV(dcdpTilde(:,1),Hp,Hps);
    dcdpy = ifiltV(dcdpTilde(:,2),Hp,Hps);
    dcdpz = ifiltV(dcdpTilde(:,3),Hp,Hps);

    f0val  = scale*compliance(iter);
    df0dxy = [dcdx(actv);dcdy];

    %% Volume constraint sensitivity
    dgvdxPhys = ones(nelz,nely,nelx);
    dgvdxPhys(Passive ~= 0) = 0;

    dgvdx = ifiltV( ...
        dgvdxPhys(:).*dx/nElem, ...
        Hx,Hxs);

    dgvdy = zeros(nElem,1);

    %% Fiber-content constraint sensitivity
    dgfdxPhysV = yPhys(:);
    dgfdyPhysV = xPhys(:);

    dgfdxPhysV(Passive(:) ~= 0) = 0;

    dgfdx = ifiltV( ...
        dgfdxPhysV.*dx/nElem, ...
        Hx,Hxs);

    dgfdy = ifiltV( ...
        dgfdyPhysV/nElem, ...
        Hy,Hys);

    fval = zeros(Nc,1);

    fval(1,1) = ...
        mean(xPhys(:))/volfrac - 1;

    fval(2,1) = ...
        mean(xPhys(:).*yPhys(:))/contentfrac - 1;

    dfvaldxy = zeros(Nc,Nd);

    dfvaldxy(1,:) = ...
        [dgvdx(actv)',dgvdy']/volfrac;

    dfvaldxy(2,:) = ...
        [dgfdx(actv)',dgfdy']/contentfrac;

    timeHis(iter,5) = toc;

    %% MMA update of density and fiber volume fraction
    tic;

    [xymma,~,~,~,~,~,~,~,~,xyL,xyU] = mmasub( ...
        Nc,Nd,iter, ...
        xyval,xymin,xymax,xyold1,xyold2, ...
        f0val,df0dxy,fval,dfvaldxy, ...
        xyL,xyU,aa0,aa,cc,dd);

    x(actv) = xymma(1:nActv,1);

    y = reshape( ...
        xymma(nActv+1:end,1), ...
        nelz,nely,nelx);

    xyold2 = xyold1;
    xyold1 = xyval;
    xyval  = [x(actv);y(:)];

    change = max(abs(xyval-xyold1));

    timeHis(iter,6) = toc;

    %% Fiber-orientation update
    tic;

    df0dp = [dcdpx;dcdpy;dcdpz];
    pval  = [px(:);py(:);pz(:)];

    pval = updateOrientationVariables( ...
        df0dp,pval,iter);

    px = reshape( ...
        pval(1:nElem), ...
        nelz,nely,nelx);

    py = reshape( ...
        pval(nElem+(1:nElem)), ...
        nelz,nely,nelx);

    pz = reshape( ...
        pval(2*nElem+(1:nElem)), ...
        nelz,nely,nelx);

    timeHis(iter,7) = toc;

    %% Print iteration information
    fprintf([ ...
        ' It.:%5i Obj.:%11.4e Com.:%11.4e ', ...
        'Vol.:%7.3f fVf.:%7.3f ch.:%7.3f ', ...
        'MGCGit.:%5i MGCGt.:%8.4f s\n'], ...
        iter, ...
        OBJHis(iter,1), ...
        compliance(iter,1), ...
        volumefrac(iter,1), ...
        contentfracHis(iter,1), ...
        change, ...
        mgcgIterHis(iter,1), ...
        mgcgTimeHis(iter,1));

    %% Update Heaviside parameter
    if heavbeta < 128 && ...
            (loopbeta >= 30 || change < tol)

        heavbeta = 2*heavbeta;
        loopbeta = 0;
        change   = 1;

        fprintf( ...
            'Parameter heavbeta increased to %g.\n', ...
            heavbeta);
    end
end

%% Trim histories
compliance     = compliance(1:iter,1);
volumefrac     = volumefrac(1:iter,1);
contentfracHis = contentfracHis(1:iter,1);
OBJHis         = OBJHis(1:iter,1);
timeHis        = timeHis(1:iter,:);

mgcgIterHis = mgcgIterHis(1:iter,1);
mgcgTimeHis = mgcgTimeHis(1:iter,1);

% For consistency with the thermo-mechanical result files.
% Because this code has only one MGCG call per optimization iteration,
% these totals are identical to the detailed histories above.
mgcgTotalIterPerOptIter = mgcgIterHis;
mgcgTotalTimePerOptIter = mgcgTimeHis;

%% Reconstruct full stiffness matrix for validation and saving
[~,~,S0] = computeConstitutiveMatrix3D( ...
    Em,vm,Ef1,Ef2,vf12,vf23,Gf12,yPhys);

D = zeros(6,6,nElem);

for k = 1:numD
    Dk = reshape(Dcomp(k,:),1,1,nElem);

    D(rowD(k),colD(k),:) = Dk;
    D(colD(k),rowD(k),:) = Dk;
end

%% Optional visualization and validation
plotOptimizedDesignMech(he,he,he, xPhys,yPhys, pxPhys,pyPhys,pzPhys, U);
plotHistory(compliance,volumefrac,contentfracHis,'Mean Compliance');
verifyOptimality3D(he,he,he,U,D,S0, xPhys,yPhys, pxPhys,pyPhys,pzPhys);

%% Save results
filename = sprintf( ...
    ['CFRCTop3D_MeanCompliance_%dx%dx%d_', ...
     'vf%.1fcf%.2f_Fm%.1f.mat'], ...
    nelx,nely,nelz, ...
    volfrac,contentfrac,FmSum);

save(filename, ...
    'compliance', ...
    'volumefrac', ...
    'contentfracHis', ...
    'OBJHis', ...
    'xPhys', ...
    'yPhys', ...
    'pxPhys', ...
    'pyPhys', ...
    'pzPhys', ...
    'U', ...
    'D', ...
    'S0', ...
    'he', ...
    'FmSum', ...
    'timeHis', ...
    'timeLabels', ...
    'mgcgIterHis', ...
    'mgcgTimeHis', ...
    'mgcgCallLabels', ...
    'mgcgTotalIterPerOptIter', ...
    'mgcgTotalTimePerOptIter');

fprintf('\nOptimization completed.\n');
fprintf('Result file: %s\n',filename);

end
