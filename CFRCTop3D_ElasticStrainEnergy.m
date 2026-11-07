function CFRCTop3D_ElasticStrainEnergy(nelx,nely,nelz,he,volfrac,contentfrac,rxmin,rymin,rpmin,cgtol,cgmax,maxiter,FmSumInput,loadBandHalfWidthInput)
% CFRCTop3D_ElasticStrainEnergy
% -------------------------------------------------------------------------
% Mechanical-load-only topology optimization of 3D continuous
% fiber-reinforced composite structures.
%
% Design variables:
%   x  : structural density
%   y  : local fiber volume fraction
%   p  : fiber orientation vector [px,py,pz]
%
% Objective:
%   Minimum elastic strain energy under the prescribed mechanical load,
%
%       Phi = 1/2 * U' * K * U = 1/2 * Fm' * U.
%
% No thermal analysis is performed:
%   - no heat-conduction equation;
%   - no temperature field;
%   - no thermal-expansion load;
%   - no thermal adjoint equation;
%   - no thermal sensitivity terms.
%
% The support, mechanical loading region, passive solid layers, filters,
% Heaviside continuation, and design-variable update scheme are retained
% from CFRCTop3D_ThermoMech_StrainEnergy_MGPCG.
%
% The passive end layers are topologically non-designable only:
% xPhys is fixed to one there, whereas y and [px,py,pz] remain design
% variables and may vary spatially.
%
% Example:
% CFRCTop3D_ElasticStrainEnergy(80,40,40,1/80,0.30,0.09,1.5,1.5,1.5,1e-5,1000,200,-1.0e6,1);
% -------------------------------------------------------------------------

%% Default inputs
if nargin < 14 || isempty(loadBandHalfWidthInput)
    loadBandHalfWidthInput = 1;
end

if nargin < 13 || isempty(FmSumInput)
    FmSumInput = -1.0e6;
end

%% Input checks
if ~isscalar(FmSumInput) || ~isfinite(FmSumInput)
    error('FmSumInput must be a finite scalar.');
end

if ~isscalar(loadBandHalfWidthInput) || ...
        ~isfinite(loadBandHalfWidthInput) || ...
        loadBandHalfWidthInput < 0
    error('loadBandHalfWidthInput must be a finite nonnegative scalar.');
end

if volfrac <= 0 || volfrac > 1
    error('volfrac must lie in (0,1].');
end

if contentfrac <= 0
    error('contentfrac must be positive.');
end

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
Em = 4.35e9;
vm = 0.37;

Ef1  = 233.13e9;
Ef2  = 23.11e9;
vf12 = 0.20;
vf23 = 0.40;
Gf12 = 8.97e9;

ymax = 0.5;

if contentfrac > ymax*volfrac
    error(['contentfrac must satisfy contentfrac <= ymax*volfrac. ', ...
           'Current values give %.6g > %.6g.'], ...
           contentfrac,ymax*volfrac);
end

%% Support and mechanical load
nElem = nelx*nely*nelz;
nNode = (nelx+1)*(nely+1)*(nelz+1);
nDof  = 3*nNode;

% Both x-end surfaces are fully clamped.
[kfs,jfs,ifs] = ndgrid( ...
    1:nelz+1, ...
    1:nely+1, ...
    [1,nelx+1]);

fixednid = sub2ind( ...
    [nelz+1,nely+1,nelx+1], ...
    kfs(:),jfs(:),ifs(:));

fixeddofs = [ ...
    3*fixednid(:)-2; ...
    3*fixednid(:)-1; ...
    3*fixednid(:)];

% Retain exactly the same mechanical-load construction as the
% thermo-mechanical double-clamped-beam code.
FmSum = FmSumInput;

[Fm,loadnid,loaddofs,loadWeights,loadInfo] = ...
    buildUniformBandLoad3D( ...
        nelx,nely,nelz,nDof, ...
        FmSum,loadBandHalfWidthInput);

U = zeros(nDof,1);

%% Passive solid layers
Passive = zeros(nelz,nely,nelx);

nPassiveX = min( ...
    max(1,round(nelx/40)), ...
    floor(nelx/2));

passiveX = unique([ ...
    1:nPassiveX, ...
    nelx-nPassiveX+1:nelx]);

% Passive = 2 fixes density to solid.
% Fiber volume fraction and fiber orientation remain design variables.
Passive(:,:,passiveX) = 2;

actv  = find(Passive == 0);
nActv = numel(actv);

passiveSolidFrac = nnz(Passive == 2)/nElem;

if volfrac <= passiveSolidFrac || volfrac > 1
    error(['volfrac must be in (%.6f,1] because passive solid elements ', ...
           'are included in the volume constraint.'], passiveSolidFrac);
end

%% Finite-element preparation
nodenrs = int32(reshape( ...
    1:nNode, ...
    1+nelz,1+nely,1+nelx));

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

% Lower-triangular sparse-assembly indices
[sI,sII] = find(tril(ones(24)));
[iK,jK] = deal(edofMat(:,sI)',edofMat(:,sII)');

%% Mechanical template stiffness matrices
numD = 21;

KET  = zeros(24,24,numD);
KETc = zeros(300,numD);

trilIndex = tril(true(24));
[rowD,colD] = find(tril(ones(6)));

for n = 1:numD
    DT = zeros(6,6);

    i = rowD(n);
    j = colD(n);

    DT(i,j) = 1;
    DT(j,i) = 1;

    KEt = H8BasicKe(he,he,he,DT);

    KET(:,:,n) = KEt;
    KETc(:,n)  = KEt(trilIndex);
end

%% Mechanical multigrid operators
PuK = cell(nl-1,1);

for l = 1:nl-1
    PuK{l,1} = prepcoarse( ...
        nelz/2^(l-1), ...
        nely/2^(l-1), ...
        nelx/2^(l-1));
end

NK = ones(nDof,1);
NK(fixeddofs) = 0;

NullK = spdiags(NK,0,nDof,nDof);

% The mechanical right-hand side is independent of the design.
F = NullK*Fm;

%% Filtering and projection
[Hx,Hxs] = PreFILTER3D(nelx,nely,nelz,rxmin);
[Hy,Hys] = PreFILTER3D(nelx,nely,nelz,rymin);
[Hp,Hps] = PreFILTER3D(nelx,nely,nelz,rpmin);

filt = @(v,H,Hs) ...
    convn(v,H,'same')./Hs;

ifiltV = @(v,H,Hs) ...
    reshape( ...
        convn( ...
            reshape(v,nelz,nely,nelx)./Hs, ...
            H,'same'), ...
        [],1);

prj = @(v,eta,beta) ...
    (tanh(beta*eta)+tanh(beta*(v-eta)))./ ...
    (tanh(beta*eta)+tanh(beta*(1-eta)));

dprj = @(v,eta,beta) ...
    beta*(1-tanh(beta*(v-eta)).^2)./ ...
    (tanh(beta*eta)+tanh(beta*(1-eta)));

%% Initial design
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

% Heaviside projection continuation:
% beta = 1,2,4,...,128.
loopbeta = 0;
heavbeta = 1;
heaveta  = 0.5;

%% MMA variables
Nc = 2;
Nd = nActv+nElem;

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

%% Histories
strainEnergy     = zeros(maxiter,1);
compliance       = zeros(maxiter,1);
volumefrac       = zeros(maxiter,1);
contentfracHis   = zeros(maxiter,1);
OBJHis           = zeros(maxiter,1);
iterationTimeHis = zeros(maxiter,1);

% Columns:
% 1 filtering/projection/normalization
% 2 mechanical material transformation
% 3 stiffness assembly
% 4 mechanical MGPCG solution
% 5 mechanical sensitivity analysis
% 6 MMA update of density and fiber volume fraction
% 7 fiber-orientation update
timeHis = zeros(maxiter,7);

timeLabels = { ...
    'FilteringProjectionNormalization', ...
    'MechanicalMaterialTransformation', ...
    'MechanicalStiffnessAssembly', ...
    'MechanicalMGPCGSolve', ...
    'MechanicalSensitivityAnalysis', ...
    'MMAUpdateDensityFiberFraction', ...
    'FiberOrientationUpdate'};

% Detailed MGCG history.
% This mechanical-only formulation calls MGCG once per optimization
% iteration to solve the mechanical equilibrium equation K*U = F.
mgcgIterHis = zeros(maxiter,1);
mgcgTimeHis = zeros(maxiter,1);
mgcgCallLabels = {'Mechanical equilibrium'};

totalTimer = tic;

%% Optimization loop
while change > tol && iter < maxiter
    iterationTimer = tic;

    iter = iter+1;
    loopbeta = loopbeta+1;

    %% Filtering, projection, and orientation normalization
    stageTimer = tic;

    xTilde = filt(x,Hx,Hxs);
    yPhys  = filt(y,Hy,Hys);

    xPhys = prj(xTilde,heaveta,heavbeta);

    xPhys(Passive == 1) = 0;
    xPhys(Passive == 2) = 1;

    pxTilde = filt(px,Hp,Hps);
    pyTilde = filt(py,Hp,Hps);
    pzTilde = filt(pz,Hp,Hps);

    Lp = sqrt( ...
        pxTilde.^2 + ...
        pyTilde.^2 + ...
        pzTilde.^2);

    Lp = max(Lp,1e-12);

    pxPhys = pxTilde./Lp;
    pyPhys = pyTilde./Lp;
    pzPhys = pzTilde./Lp;

    timeHis(iter,1) = toc(stageTimer);

    %% Mechanical material properties
    stageTimer = tic;

    l1 = pxPhys(:);
    m1 = pyPhys(:);
    n1 = pzPhys(:);

    [xmpD,dxmpD] = MaterialInterpolationD(xPhys);

    [D0,dD0y,~] = computeConstitutiveMatrix3D( ...
        Em,vm, ...
        Ef1,Ef2,vf12,vf23,Gf12, ...
        yPhys);

    [Dcomp,dDlcomp,dDmcomp,dDncomp,dDycomp] = ...
        ComputeElasticityMatrices3DLowerTri21( ...
            D0,l1,m1,n1,dD0y);

    timeHis(iter,2) = toc(stageTimer);

    %% Assemble mechanical stiffness entries
    stageTimer = tic;

    sK = KETc*(Dcomp.*xmpD(:).');

    timeHis(iter,3) = toc(stageTimer);

    %% Solve mechanical equilibrium
    stageTimer = tic;

    K = cell(nl,1);

    K{1,1} = sparse( ...
        iK,jK,sK(:), ...
        nDof,nDof);

    K{1,1} = ...
        K{1,1}+K{1,1}'-diag(diag(K{1,1}));

    clear sK

    K{1,1} = ...
        NullK'*K{1,1}*NullK - ...
        (NullK-speye(nDof,nDof));

    for l = 1:nl-1
        K{l+1,1} = ...
            PuK{l,1}'*(K{l,1}*PuK{l,1});
    end

    Lfac = chol(K{nl,1},'lower');
    Ufac = Lfac';

    mgcgTimer = tic;

    [mgcgIterHis(iter,1),~,U] = mgcg( ...
        K,F,U, ...
        Lfac,Ufac,PuK,nl,1, ...
        cgtol,cgmax);

    % Pure MGCG solution time. This excludes stiffness assembly,
    % multigrid hierarchy construction, and coarse-grid factorization.
    mgcgTimeHis(iter,1) = toc(mgcgTimer);

    % Full mechanical-equilibrium stage time.
    timeHis(iter,4) = toc(stageTimer);

    %% Objective and constraints
    compliance(iter,1) = F'*U;

    strainEnergy(iter,1) = ...
        0.5*compliance(iter,1);

    Phi = strainEnergy(iter,1);

    OBJHis(iter,1) = Phi;

    % The prescribed solid regions are included in the total volume budget.
    volumefrac(iter,1) = mean(xPhys(:));

    % The fiber-content constraint includes all solid and design elements.
    % y remains designable in the passive solid end layers.
    contentfracHis(iter,1) = ...
        mean(xPhys(:).*yPhys(:));

    %% Mechanical sensitivity analysis
    stageTimer = tic;

    if mod(iter,20) == 1
        scale = 1/max(abs(Phi),eps);
    end

    Ud = U(edofMat);

    % uku(:,n) = -u_e' * KET_n * u_e
    uku = zeros(nElem,numD);

    for n = 1:numD
        uku(:,n) = ...
            -sum((Ud*KET(:,:,n)).*Ud,2);
    end

    % Phi = 1/2 F' U:
    % dPhi/dz = -1/2 U' (dK/dz) U.
    cex = 0.5*dxmpD(:).* ...
        sum(Dcomp.'.*uku,2);

    cey = 0.5*xmpD(:).* ...
        sum(dDycomp.'.*uku,2);

    cepx = 0.5*xmpD(:).* ...
        sum(dDlcomp.'.*uku,2);

    cepy = 0.5*xmpD(:).* ...
        sum(dDmcomp.'.*uku,2);

    cepz = 0.5*xmpD(:).* ...
        sum(dDncomp.'.*uku,2);

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

    % Density is fixed in passive regions.
    dcdxPhysV(Passive(:) ~= 0) = 0;

    %% Filter objective sensitivities
    dx = dprj( ...
        xTilde(:), ...
        heaveta,heavbeta);

    dcdx = ifiltV( ...
        dcdxPhysV.*dx, ...
        Hx,Hxs);

    dcdy = ifiltV( ...
        dcdyPhysV, ...
        Hy,Hys);

    % Chain rule through orientation normalization.
    G = [ ...
        dcdpxPhysV, ...
        dcdpyPhysV, ...
        dcdpzPhysV];

    P = [l1,m1,n1];

    dcdpTilde = ...
        (G-P.*sum(P.*G,2))./Lp(:);

    dcdpx = ifiltV( ...
        dcdpTilde(:,1), ...
        Hp,Hps);

    dcdpy = ifiltV( ...
        dcdpTilde(:,2), ...
        Hp,Hps);

    dcdpz = ifiltV( ...
        dcdpTilde(:,3), ...
        Hp,Hps);

    f0val  = scale*Phi;
    df0dxy = [dcdx(actv);dcdy(:)];

    %% Constraint sensitivities
    % Total structural-volume constraint; passive solid regions count
    % toward the volume budget but have zero density sensitivity.
    dgvdxPhys = ones(nelz,nely,nelx);
    dgvdxPhys(Passive ~= 0) = 0;

    dgvdx = ifiltV( ...
        dgvdxPhys(:).*dx/nElem, ...
        Hx,Hxs);

    dgvdy = zeros(nElem,1);

    % Global fiber-content constraint
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
        mean(xPhys(:))/volfrac-1;

    fval(2,1) = ...
        mean(xPhys(:).*yPhys(:))/contentfrac-1;

    dfvaldxy = zeros(Nc,Nd);

    dfvaldxy(1,:) = ...
        [dgvdx(actv)',dgvdy(:)']/volfrac;

    dfvaldxy(2,:) = ...
        [dgfdx(actv)',dgfdy(:)']/contentfrac;

    timeHis(iter,5) = toc(stageTimer);

    %% MMA update of x and y
    stageTimer = tic;

    [xymma,~,~,~,~,~,~,~,~,xyL,xyU] = ...
        mmasub( ...
            Nc,Nd,iter, ...
            xyval,xymin,xymax, ...
            xyold1,xyold2, ...
            f0val,df0dxy, ...
            fval,dfvaldxy, ...
            xyL,xyU, ...
            aa0,aa,cc,dd);

    x(actv) = xymma(1:nActv,1);

    y = reshape( ...
        xymma(nActv+1:end,1), ...
        nelz,nely,nelx);

    xyold2 = xyold1;
    xyold1 = xyval;
    xyval  = [x(actv);y(:)];

    change = max(abs(xyval-xyold1));

    timeHis(iter,6) = toc(stageTimer);

    %% Update fiber orientations
    stageTimer = tic;

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

    timeHis(iter,7) = toc(stageTimer);

    iterationTimeHis(iter,1) = ...
        toc(iterationTimer);

    %% Print results
    fprintf([ ...
        ' It.:%5i Obj.:%11.4e Phi:%11.4e ', ...
        'Com.:%11.4e Vol.:%7.3f fVf.:%7.3f ', ...
        'ch.:%7.3f Time:%9.3f s ', ...
        'MGCGit.:%5i MGCGt.:%8.4f s\n'], ...
        iter, ...
        OBJHis(iter,1), ...
        strainEnergy(iter,1), ...
        compliance(iter,1), ...
        volumefrac(iter,1), ...
        contentfracHis(iter,1), ...
        change, ...
        iterationTimeHis(iter,1), ...
        mgcgIterHis(iter,1), ...
        mgcgTimeHis(iter,1));

    %% Heaviside continuation
    if heavbeta < 128 && ...
            (loopbeta >= 30 || change < tol)

        heavbeta = 2*heavbeta;
        loopbeta = 0;
        change   = 1;

        fprintf( ...
            'Projection parameter beta increased to %g.\n', ...
            heavbeta);
    end
end

totalWallTime = toc(totalTimer);

%% Trim histories
strainEnergy     = strainEnergy(1:iter,1);
compliance       = compliance(1:iter,1);
volumefrac       = volumefrac(1:iter,1);
contentfracHis   = contentfracHis(1:iter,1);
OBJHis           = OBJHis(1:iter,1);
timeHis          = timeHis(1:iter,:);
iterationTimeHis = iterationTimeHis(1:iter,1);

mgcgIterHis = mgcgIterHis(1:iter,1);
mgcgTimeHis = mgcgTimeHis(1:iter,1);

% Retain the same variable naming convention as the thermo-mechanical
% programs. Because there is only one MGCG call in each optimization
% iteration, the totals equal the detailed histories.
mgcgTotalIterPerOptIter = mgcgIterHis;
mgcgTotalTimePerOptIter = mgcgTimeHis;

averageIterationTime = mean(iterationTimeHis);

%% Reconstruct full stiffness matrix for validation
[~,~,S0] = computeConstitutiveMatrix3D( ...
    Em,vm, ...
    Ef1,Ef2,vf12,vf23,Gf12, ...
    yPhys);

D = zeros(6,6,nElem);

for k = 1:numD
    Dk = reshape(Dcomp(k,:),1,1,nElem);

    D(rowD(k),colD(k),:) = Dk;
    D(colD(k),rowD(k),:) = Dk;
end

%% Optional visualization and validation
% plotOptimizedDesignMech( ...
%     he,he,he, ...
%     xPhys,yPhys, ...
%     pxPhys,pyPhys,pzPhys, ...
%     U);
%
% plotHistory( ...
%     strainEnergy, ...
%     volumefrac, ...
%     contentfracHis, ...
%     'Elastic Strain Energy');
%
% verifyOptimality3D( ...
%     he,he,he,U,D,S0, ...
%     xPhys,yPhys, ...
%     pxPhys,pyPhys,pzPhys);

%% Save results
filename = sprintf( ...
    ['CFRCTop3D_ElasticStrainEnergy_', ...
     '%dx%dx%d_vf%.1fcf%.2f_Fm%.1f.mat'], ...
    nelx,nely,nelz, ...
    volfrac,contentfrac,FmSum);

save(filename, ...
    'strainEnergy', ...
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
    'F', ...
    'D', ...
    'S0', ...
    'he', ...
    'FmSum', ...
    'loadBandHalfWidthInput', ...
    'loadnid', ...
    'loaddofs', ...
    'loadWeights', ...
    'loadInfo', ...
    'Passive', ...
    'actv', ...
    'nPassiveX', ...
    'passiveX', ...
    'timeHis', ...
    'timeLabels', ...
    'mgcgIterHis', ...
    'mgcgTimeHis', ...
    'mgcgCallLabels', ...
    'mgcgTotalIterPerOptIter', ...
    'mgcgTotalTimePerOptIter', ...
    'iterationTimeHis', ...
    'totalWallTime', ...
    'averageIterationTime', ...
    'heavbeta', ...
    'heaveta', ...
    'rxmin', ...
    'rymin', ...
    'rpmin');

fprintf('\nMechanical-only optimization completed.\n');
fprintf('Iterations              : %d\n',iter);
fprintf('Final strain energy     : %.6e J\n',strainEnergy(end));
fprintf('Final compliance        : %.6e J\n',compliance(end));
fprintf('Total wall-clock time   : %.3f s\n',totalWallTime);
fprintf('Average time/iteration  : %.3f s\n',averageIterationTime);
fprintf('Total MGCG iterations   : %d\n',sum(mgcgIterHis));
fprintf('Total MGCG solve time   : %.3f s\n',sum(mgcgTimeHis));
fprintf('Average MGCG iterations : %.3f\n',mean(mgcgIterHis));
fprintf('Average MGCG solve time : %.6f s\n',mean(mgcgTimeHis));
fprintf('Result file             : %s\n',filename);

end
