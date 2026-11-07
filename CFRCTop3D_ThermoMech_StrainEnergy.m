function CFRCTop3D_ThermoMech_StrainEnergy(nelx,nely,nelz,he,volfrac,contentfrac,rxmin,rymin,rpmin,cgtol,cgmax,maxiter,qInput,FmSumInput,loadBandHalfWidthInput)
% Decoupled update method: [x;y] are updated first, then [px;py;pz] are updated separately.
% Thermo-mechanical coupled topology optimization of continuous fiber reinforced composites.
% Objective: Elastic strain energy (complementary energy form).
% Records the iteration count, final relative residual, and wall-clock
% time of every MGCG call during each optimization iteration.
% qInput is the prescribed surface heat flux on the z+ face.
% FmSumInput is the total mechanical load applied to loaddofs.
% loadBandHalfWidthInput controls the x-direction half-width of the load band.
%
% Phi = 1/2 * U^T * K_m * U - U^T * K_uT * DeltaT + Phi^ith
%     = 1/2 * (U^m)^T * F_m - 1/2 * (U^th)^T * F_th + Phi^ith
%
% where U^m and U^th are displacement components from mechanical and thermal loads.
if nargin < 15 || isempty(loadBandHalfWidthInput), loadBandHalfWidthInput = 1; end
if nargin < 14 || isempty(FmSumInput), FmSumInput = -1000000; end
if nargin < 13 || isempty(qInput), qInput = 5.0; end
if ~isscalar(loadBandHalfWidthInput), error('loadBandHalfWidthInput must be a scalar.'); end
if ~isscalar(FmSumInput), error('FmSumInput must be a scalar total mechanical load.'); end
if ~isscalar(qInput), error('qInput must be a scalar heat flux applied on the z+ face.'); end

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

%% Define support and load
nElem = nelx*nely*nelz; % Number of elements
nNode = (nelx+1)*(nely+1)*(nelz+1); % number of nodes
nDof = 3*nNode; % Number of dofs
% USER-DEFINED SUPPORT FIXED DOFs
[kfs,jfs,ifs] = ndgrid(1:nelz+1,1:nely+1,[1,nelx+1]);% Coordinates
fixednid = sub2ind([nelz+1,nely+1,nelx+1],kfs(:),jfs(:),ifs(:));% Fixed Node IDs
fixeddofs = [3*fixednid(:)-2;3*fixednid(:)-1;3*fixednid(:)]; % Fixed DOFs
% USER-DEFINED LOAD DOFs
FmSum = FmSumInput;
[Fm, loadnid, loaddofs, loadWeights, loadInfo] = buildUniformBandLoad3D(nelx,nely,nelz,nDof,FmSum,loadBandHalfWidthInput);
Um = zeros(nDof,1);  % mechanical displacement
Uth = zeros(nDof,1); % thermal displacement
U = zeros(nDof,1);
PTOri = zeros(nNode,1);
lambdaT = zeros(nNode,1); % thermal adjoint for strain energy

%% Define heat flux and thermal boundary
% 6-component vector [q1, q2, q3, q4, q5, q6] (W/m^2), corresponding to:
% 1:x-, 2:x+, 3:y-, 4:y+, 5:z-, 6:z+ face
q = qInput;
q1 = [0,0,0,0,0,q];
[kq1,jq1,iq1] = ndgrid(nelz,1:nely,1:nelx); % Elements subject to heat flux
qelem1 = sub2ind([nelz,nely,nelx],kq1(:),jq1(:),iq1(:));

BoundT = 20; % Prescribed/initial temperature on thermal boundary nodes
Tini = 20;   % Reference temperature for thermal strain
[kBN,jBN,iBN] = ndgrid(1,1:nely+1,1:nelx+1);
BoundNode = sub2ind([nelz+1,nely+1,nelx+1],kBN(:),jBN(:),iBN(:));

%% Passive elements definition
Passive = zeros(nelz,nely,nelx);
nPassiveX = min(max(1,round(nelx/40)),floor(nelx/2));
passiveX = unique([1:nPassiveX, nelx-nPassiveX+1:nelx]);
Passive(:,:,passiveX) = 2; % 1 for void, 2 for solid (non-designable domain)
% Extract active design domain indices
actv = find(Passive == 0); % Find linear indices of all designable elements
nActv = length(actv);      % Number of active variables

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
[sIM, sIIM] = find(ones(8,8));
[iM,jM] = deal(nids(:,sIM)',nids(:,sIIM)');

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

% Precompute thermal mass matrix M_T = integral N^T N dOmega (same for all cubic elements)
% Used for evaluating eta_e = DeltaT_e^T * M_T * DeltaT_e
MT = computeThermalMassMatrix(he, he, he);

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
if xinitActv < 0 || xinitActv > 1
    error('Prescribed total structural volume fraction is incompatible with the passive solid volume.');
end
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

% MMA variables: density only in the design domain, fiber volume fraction in all elements
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

strainEnergy = zeros(maxiter,1);
volumefrac = zeros(maxiter,1);
contentfracHis = zeros(maxiter,1);
OBJHis = zeros(maxiter,1);
timeHis = zeros(maxiter,8);

% MGCG call histories. Columns correspond to:
%   1) thermal state equation:       K_th * T       = P_T
%   2) mechanical state equation:    K_m  * U^m     = F_m
%   3) thermal-load displacement:    K_m  * U^th    = F_th
%   4) thermal adjoint equation:      K_th * lambda = b_adj
mgcgIterHis   = zeros(maxiter,4);
mgcgTimeHis   = zeros(maxiter,4);
mgcgRelresHis = nan(maxiter,4);

mgcgCallLabels = { ...
    'ThermalState', ...
    'MechanicalState', ...
    'ThermalDisplacement', ...
    'ThermalAdjoint'};

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

    % Compute alpha^T * D * alpha (phi0) for initial thermal energy
    % phi0e = alpha_e^T * D_e * alpha_e  -> scalar per element
    phi0 = computeQuadraticFormAlphaDalpha(Dcomp, Acomp);
    dphi0y = computeQuadraticFormAlphaDalpha(dDycomp, Acomp); % alpha^T * dDycomp * alpha
    dphi0l = computeQuadraticFormAlphaDalpha(dDlcomp, Acomp);
    dphi0m = computeQuadraticFormAlphaDalpha(dDmcomp, Acomp);
    dphi0n = computeQuadraticFormAlphaDalpha(dDncomp, Acomp);

    BetaDotAlpha = sum(Beta .* Acomp, 1)'; % nElem x 1: Beta^T * alpha
    dBetayDotAlpha = sum(dBetay .* Acomp, 1)'; % nElem x 1
    dBetalDotAlpha = sum(dBetal .* Acomp, 1)'; % nElem x 1
    dBetamDotAlpha = sum(dBetam .* Acomp, 1)';
    dBetanDotAlpha = sum(dBetan .* Acomp, 1)';

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

    mgcgTimer = tic;
    [mgcgIterHis(iter,1), mgcgRelresHis(iter,1), T] = ...
        mgcg(Kth, PT, T, Lfacth, Ufacth, PuKth, nl, 1, cgtol, cgmax);
    mgcgTimeHis(iter,1) = toc(mgcgTimer);

    timeHis(iter,4) = toc(thermalStageTimer);

    %% SOLVE MECHANICAL SYSTEM: decompose into U^m and U^th
    mechanicalStageTimer = tic;
    deltaT = T - Tini;
    deltaTe = deltaT(nids)'; % 8 x nElem matrix

    % Assemble thermal load vector Fth
    scaledBeta = Beta .* xmpBeta(:).';
    V3D = reshape(deltaTe, 8, 1, []) .* reshape(scaledBeta, 1, 6, []);
    sFth = HET2D * reshape(V3D, 48, []);
    sH = HETc * (Beta .* xmpBeta(:).');
    Fth = sparse(iFth,1,sFth(:),nDof,1);

    % Assemble stiffness matrix
    K = cell(nl,1);
    K{1,1} = sparse(iK,jK,sK(:),nDof,nDof);
    K{1,1} = K{1,1}+K{1,1}'-diag(diag(K{1,1}));
    K{1,1} = NullK'*K{1,1}*NullK - (NullK-speye(nDof,nDof));

    for l = 1:nl-1
        K{l+1,1} = PuK{l,1}'*(K{l,1}*PuK{l,1});
    end
    Lfac = chol(K{nl,1},'lower'); Ufac = Lfac';

    % Solve for mechanical displacement U^m: K_m * U^m = F_m
    FmNull = NullK * Fm;

    mgcgTimer = tic;
    [mgcgIterHis(iter,2), mgcgRelresHis(iter,2), Um] = ...
        mgcg(K, FmNull, Um, Lfac, Ufac, PuK, nl, 1, cgtol, cgmax);
    mgcgTimeHis(iter,2) = toc(mgcgTimer);

    % Solve for thermal displacement U^th: K_m * U^th = F_th
    FthNull = NullK * Fth;

    mgcgTimer = tic;
    [mgcgIterHis(iter,3), mgcgRelresHis(iter,3), Uth] = ...
        mgcg(K, FthNull, Uth, Lfac, Ufac, PuK, nl, 1, cgtol, cgmax);
    mgcgTimeHis(iter,3) = toc(mgcgTimer);

    % Total displacement (for visualization)
    U = Um + Uth;
    F = Fm + Fth;

    timeHis(iter,5) = toc(mechanicalStageTimer);

    %% COMPUTE STRAIN ENERGY OBJECTIVE
    PhiMech = 0.5 * U' * F;
    PhiTherm = U' * Fth;

    % Compute initial thermal strain energy Phi^ith
    % Phi^ith = sum_e 0.5 * eta_e * alpha_e^T * D_e * alpha_e
    etae = computeEtaPerElement(MT, deltaTe);
    % Penalty factor: alpha = (x_beta/x_D)*alpha_0, D = x_D*D_0
    % => alpha^T*D*alpha = (x_beta^2/x_D) * alpha_0^T*D_0*alpha_0
    penalIth = (xmpBeta(:).^2) ./ xmpD(:);
    PhiIth = 0.5 * sum(penalIth .* etae(:) .* phi0(:));

    Phi = PhiMech - PhiTherm + PhiIth;

    strainEnergy(iter,1) = Phi;
    OBJHis(iter,1) = Phi;
    volumefrac(iter,1) = mean(xPhys(:));
    contentfracHis(iter,1) = mean(xPhys(:).*yPhys(:));

    %% SENSITIVITY ANALYSIS
    sensitivityStageTimer = tic;
    % --- Step 1: Compute thermal adjoint lambda_T^Phi ---
    % Solve: K_th * lambda_T = K_uT^T * U^th - dPhi^ith/dT
    % Element contribution from K_uT^T * U^th: bKuT_e = H_e^T * u_e^th
    Hth = sparse(double(iH), double(jH), sH(:), nDof, nNode);
    bKuT = Hth' * Uth;

    % Element contribution from dPhi^ith/dT: bPhi_e = Penalphi0_e * M_T * deltaT_e
    % The temperature derivative of the initial thermal strain energy is also scaled by the material interpolation.
    phi0penalized = penalIth .* phi0(:);
    sM = MT(:) * phi0penalized'; % Size: 64 x nElem
    Mthphi = sparse(double(iM), double(jM), sM(:), nNode, nNode);
    bPhi = Mthphi * deltaT; % deltaT = T - Tini has already been evaluated above.

    % Assemble thermal adjoint RHS b = bKuT - bPhi
    bAdj = bKuT - bPhi;
    bAdj(BoundNode) = 0; % Apply the thermal adjoint boundary conditions.

    % Solve for thermal adjoint
    mgcgTimer = tic;
    [mgcgIterHis(iter,4), mgcgRelresHis(iter,4), lambdaT] = ...
        mgcg(Kth, bAdj, lambdaT, Lfacth, Ufacth, PuKth, nl, 1, cgtol, cgmax);
    mgcgTimeHis(iter,4) = toc(mgcgTimer);

    %% --- Step 2: Compute element-level sensitivity contributions ---
    if mod(iter,20) == 1
        scale = 1/max(abs(Phi),eps);
    end

    % Extract element displacement vectors
    Udm  = Um(edofMat);  % nElem x 24
    Udth = Uth(edofMat); % nElem x 24

    % Compute uku (u^T * Ke * u) for both displacement fields
    ukum  = zeros(nElem, numD);
    uthkuth = zeros(nElem, numD);
    for n = 1:numD
        ukum(:,n)  = -sum((Udm  * KET(:,:,n)) .* Udm,  2);
        uthkuth(:,n) = -sum((Udth * KET(:,:,n)) .* Udth, 2);
    end

    % Compute uFth_th for coupling sensitivity
    % uthFth(e,n) = 2 * (u_e^th)^T * HET(:,:,n) * deltaT_e
    uthFth = zeros(nElem, numH);
    for n = 1:numH
        uthFth(:,n) = 2 * sum((Udth * HET(:,:,n)) .* deltaTe', 2);
    end

    % Compute thermal adjoint contribution: akthT_phi(e,n) = lambda_T_e^T * KthET_n * T_e
    lambdaTN = lambdaT(nids);  % nElem x 8
    TN = T(nids);             % nElem x 8
    akthTphi = zeros(nElem, numKa);
    for n = 1:numKa
        akthTphi(:,n) = sum((lambdaTN * KthET(:,:,n)) .* TN, 2);
    end

    %% --- Step 3: Assemble x-sensitivity (density) ---
    % Terms 1+3: mechanical stiffness contributions with U^m and U^th
    cex = 0.5 * dxmpD(:) .* (sum(Dcomp.' .* ukum, 2) - sum(Dcomp.' .* uthkuth, 2));

    % Term 2: thermo-mechanical coupling with U^th
    cex = cex - 0.5 * dxmpBeta(:) .* sum(Beta.' .* uthFth, 2);

    % Term 4: initial thermal energy derivative w.r.t. x
    % dPhi^ith/dx = eta_e * dxmpBeta * Beta^T*alpha - 0.5 * eta_e * dxmpD * (alpha^T*Dcomp*alpha)
    
    % Apply the chain rule consistently for alpha_e = (x_beta / x_D) * alpha_0.
    factor1_x = dxmpBeta(:) .* (xmpBeta(:) ./ xmpD(:));
    factor2_x = 0.5 * dxmpD(:) .* ((xmpBeta(:) ./ xmpD(:)).^2);
    cex = cex + factor1_x .* etae(:) .* BetaDotAlpha(:);
    cex = cex - factor2_x .* etae(:) .* phi0(:);

    % Term 5: thermal adjoint contribution lambda_T^T * dKth/dx * T
    cex = cex + dxmpKa(:) .* sum(Kacomp.' .* akthTphi, 2);

    cex = scale * cex;

    %% --- Step 4: Assemble y-sensitivity (fiber volume fraction) ---
    % Terms 1+3: mechanical stiffness with dD/dy
    cey = 0.5 * xmpD(:) .* (sum(dDycomp.' .* ukum, 2) - sum(dDycomp.' .* uthkuth, 2));

    % Term 2: thermo-mechanical coupling with dBeta/dy
    cey = cey - 0.5 * xmpBeta(:) .* sum(dBetay.' .* uthFth, 2);

    % Term 4: initial thermal energy derivative w.r.t. y
    % dPhi^ith/dy = eta_e * dBetay^T*alpha - 0.5 * eta_e * alpha^T*dDycomp*alpha
    
    % Apply the common interpolation scaling factor penalIth to the initial thermal energy terms.
    cey = cey + penalIth .* etae(:) .* dBetayDotAlpha(:);
    cey = cey - 0.5 * penalIth .* etae(:) .* dphi0y(:);

    % Term 5: thermal adjoint with dKth/dy
    cey = cey + xmpKa(:) .* sum(dKaycomp.' .* akthTphi, 2);

    cey = scale * cey;

    %% --- Step 5: Assemble orientation sensitivity (l, m, n) ---
    % Directional derivatives of strain energy

    % Terms 1+3: mechanical stiffness directional derivatives
    cepx = 0.5 * xmpD(:) .* (sum(dDlcomp.' .* ukum, 2) - sum(dDlcomp.' .* uthkuth, 2));
    cepy = 0.5 * xmpD(:) .* (sum(dDmcomp.' .* ukum, 2) - sum(dDmcomp.' .* uthkuth, 2));
    cepz = 0.5 * xmpD(:) .* (sum(dDncomp.' .* ukum, 2) - sum(dDncomp.' .* uthkuth, 2));

    % Term 2: thermo-mechanical coupling directional derivatives
    cepx = cepx - 0.5 * xmpBeta(:) .* sum(dBetal.' .* uthFth, 2);
    cepy = cepy - 0.5 * xmpBeta(:) .* sum(dBetam.' .* uthFth, 2);
    cepz = cepz - 0.5 * xmpBeta(:) .* sum(dBetan.' .* uthFth, 2);

    % Term 4: initial thermal energy directional derivatives
    % dPhi^ith/d(l) = eta_e * dBetal^T*alpha - 0.5*eta_e * alpha^T*dDlcomp*alpha
    

    % Apply the same penalIth interpolation scaling to the orientation derivatives.
    cepx = cepx + penalIth .* etae(:) .* dBetalDotAlpha(:) - 0.5 * penalIth .* etae(:) .* dphi0l(:);
    cepy = cepy + penalIth .* etae(:) .* dBetamDotAlpha(:) - 0.5 * penalIth .* etae(:) .* dphi0m(:);
    cepz = cepz + penalIth .* etae(:) .* dBetanDotAlpha(:) - 0.5 * penalIth .* etae(:) .* dphi0n(:);

    % Term 5: thermal adjoint directional derivatives
    cepx = cepx + xmpKa(:) .* sum(dKalcomp.' .* akthTphi, 2);
    cepy = cepy + xmpKa(:) .* sum(dKamcomp.' .* akthTphi, 2);
    cepz = cepz + xmpKa(:) .* sum(dKancomp.' .* akthTphi, 2);

    cepx = scale * cepx;
    cepy = scale * cepy;
    cepz = scale * cepz;

    dcdxPhysV = cex; dcdyPhysV = cey;
    dcdpxPhysV = cepx; dcdpyPhysV = cepy; dcdpzPhysV = cepz;

    %% FILTERING OF SENSITIVITIES
    % Objective function
    dx = dprj(xTilde(:), heaveta, heavbeta);
    dcdxPhysV(Passive(:)~=0) = 0;
    dcdx = ifiltV(dcdxPhysV .* dx, Hx, Hxs);
    dcdy = ifiltV(dcdyPhysV, Hy, Hys);

    G = [dcdpxPhysV, dcdpyPhysV, dcdpzPhysV];
    P = [l1,m1,n1];
    dcdpTiles = (G-P.*sum(P.*G,2))./Lp(:);
    dcdpx = ifiltV(dcdpTiles(:,1), Hp, Hps);
    dcdpy = ifiltV(dcdpTiles(:,2), Hp, Hps);
    dcdpz = ifiltV(dcdpTiles(:,3), Hp, Hps);

    f0val = scale * Phi;
    df0dxy = [dcdx(actv); dcdy(:)];

    % Constraint functions
    dgvdxPhys = ones(nelz,nely,nelx);
    dgvdxPhys(Passive~=0) = 0;
    dgvdx = ifiltV(dgvdxPhys(:).*dx/nElem, Hx, Hxs);
    dgvdy = zeros(nElem,1);

    dgfdxPhysV = yPhys(:); dgfdyPhysV = xPhys(:);
    dgfdxPhysV(Passive(:)~=0) = 0;
    dgfdx = ifiltV(dgfdxPhysV.*dx/nElem, Hx, Hxs);
    dgfdy = ifiltV(dgfdyPhysV/nElem, Hy, Hys);

    fval = zeros(Nc,1);
    fval(1,1) = mean(xPhys(:))/volfrac - 1;
    fval(2,1) = mean(xPhys(:).*yPhys(:))/contentfrac - 1;

    dfvaldxy = zeros(Nc,Nd);
    dfvaldxy(1,:) = [dgvdx(actv)', dgvdy(:)'] / volfrac;
    dfvaldxy(2,:) = [dgfdx(actv)', dgfdy(:)'] / contentfrac;

    timeHis(iter,6) = toc(sensitivityStageTimer);

    %% DESIGN VARIABLE UPDATE
    % Update active density variables and global fiber volume fraction variables
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
    fprintf([' It.:%5i Obj.:%11.4f Phi:%11.4f Vol.:%7.3f ', ...
             'fVf.:%7.3f ch.:%7.3f MGCGit.:%5d MGCGt.:%9.4f s\n'], ...
        iter, OBJHis(iter,1), strainEnergy(iter,1), volumefrac(iter,1), ...
        contentfracHis(iter,1), change, ...
        sum(mgcgIterHis(iter,:)), sum(mgcgTimeHis(iter,:)));

    %% UPDATE HEAVISIDE REGULARIZATION PARAMETER
    if heavbeta < 128 && (loopbeta >= 30 || change < tol)
        heavbeta = 2*heavbeta;
        loopbeta = 0;
        change = 1;
        fprintf('Parameter heavbeta increased to %g.\n', heavbeta);
    end
end

strainEnergy = strainEnergy(1:iter,1);
volumefrac = volumefrac(1:iter,1);
contentfracHis = contentfracHis(1:iter,1);
timeHis = timeHis(1:iter,:);

mgcgIterHis   = mgcgIterHis(1:iter,:);
mgcgTimeHis   = mgcgTimeHis(1:iter,:);
mgcgRelresHis = mgcgRelresHis(1:iter,:);

mgcgTotalIterPerOptIter = sum(mgcgIterHis,2);
mgcgTotalTimePerOptIter = sum(mgcgTimeHis,2);

%% VISUALIZATION OF THE OPTIMIZED DESIGN
% plotOptimizedDesign(he,he,he,xPhys,yPhys,pxPhys,pyPhys,pzPhys,U,T)
% plotHistory(strainEnergy,volumefrac,contentfracHis,'Elastic Strain Energy')

%% VALIDATION OF THE OPTIMALITY OF FIBER ORIENTATION
[~,~,S0] = compositeMixingParams(km1,km2,am,Em,vm,kf1,kf2,af1,af2,Ef1,Ef2,vf12,vf23,Gf12,yPhys);
% Reconstruct the full 6-by-6-by-nElem stiffness matrix
D = zeros(6,6,nElem);
for k = 1:numD
    Dk = reshape(Dcomp(k,:),1,1,nElem);
    D(rowD(k),colD(k),:) = Dk;
    D(colD(k),rowD(k),:) = Dk;
end
% verifyOptimality3D(he,he,he,U,D,S0,xPhys,yPhys,pxPhys,pyPhys,pzPhys);
filename = sprintf(['CFRCTop3D_ThermoMech_StrainEnergy_MGPCG_MGCGHistory_', ...
    '%dx%dx%d_vf%.1fcf%.2f_q%.1fFm%.1f.mat'], ...
    nelx,nely,nelz,volfrac,contentfrac,q,FmSum);
save(filename, ...
    'strainEnergy','volumefrac','contentfracHis', ...
    'xPhys','yPhys','pxPhys','pyPhys','pzPhys', ...
    'U','D','S0','he','T','q','FmSum', ...
    'loadnid','loaddofs','loadWeights','loadInfo', ...
    'Passive','actv','nPassiveX','timeHis', ...
    'mgcgIterHis','mgcgTimeHis','mgcgRelresHis', ...
    'mgcgCallLabels', ...
    'mgcgTotalIterPerOptIter','mgcgTotalTimePerOptIter');

fprintf('\nMGCG statistics saved to: %s\n',filename);
fprintf('Total MGCG calls          : %d\n',numel(mgcgIterHis));
fprintf('Total MGCG iterations     : %d\n',sum(mgcgIterHis,'all'));
fprintf('Total MGCG solve time     : %.6f s\n',sum(mgcgTimeHis,'all'));

% =========================================================================
% Nested helper functions
% =========================================================================

end
