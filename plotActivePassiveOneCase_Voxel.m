function plotActivePassiveOneCase_Voxel(matFile)
% plotActivePassiveOneCase_Voxel
% Single-case interactive 3D visualization using the ORIGINAL voxel geometry.
% Active region: opaque voxel surface with thin edge lines.
% Passive region: translucent box (top passive solid layers).
% Only a MATLAB .fig is saved; no PNG is generated.
%
% Usage:
%   plotActivePassiveOneCase_Voxel
%   plotActivePassiveOneCase_Voxel('D:\results\...q0.0Fm-1000000.0.mat')

baseDir = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(matFile)
    candidates = dir(fullfile(baseDir, 'CFRCTop3D_ThermoMech_MGPCG_MGCGHistory_*_q2.0Fm*.mat'));
    if isempty(candidates)
        [name, folder] = uigetfile('*.mat', 'Select ONE topology result (.mat)');
        if isequal(name,0), return; end
        matFile = fullfile(folder, name);
    else
        matFile = fullfile(candidates(1).folder, candidates(1).name);
    end
end
if ~isfile(matFile)
    error('MAT file not found: %s', matFile);
end

% ------------------- editable settings -------------------
opt.threshold     = 0.50;
opt.activeColor   = [0.78, 0.80, 0.82];
opt.edgeColor     = [0.58, 0.60, 0.63];
opt.edgeAlpha     = 0.38;
opt.edgeWidth     = 0.18;
opt.passiveColor  = [0.88, 0.89, 0.90];
opt.passiveAlpha  = 0.14;
opt.background    = [0.965, 0.965, 0.970];
opt.azimuth       = 30;
opt.elevation     = 30;     % keep consistent with existing paper figures
opt.margin        = 0.07;   % fraction of full domain on each side

% Lighting / material (voxel rendering)
opt.ambient       = 0.28;
opt.diffuse       = 0.72;
opt.specular      = 0.03;
opt.specPower     = 10;
% ---------------------------------------------------------

% Load data.
S = load(matFile,'xPhys','he');
info = whos('-file',matFile);
if any(strcmp({info.name},'Passive'))
    tmp = load(matFile,'Passive');
    S.Passive = tmp.Passive;
end
if ~isfield(S,'xPhys') || ~isfield(S,'he')
    error('Input MAT file must contain xPhys and he.');
end
rho = double(S.xPhys);
he  = double(S.he);
[nz,ny,nx] = size(rho);

% Passive region.
if isfield(S,'Passive')
    P = S.Passive;
else
    P = zeros(nz,ny,nx,'uint8');
    thickness = max(1, round(nz/20));
    P(end-thickness+1:end,:,:) = 2;
    fprintf('Passive reconstructed as top %d layers.\n', thickness);
end
assert(isequal(size(P),size(rho)), 'Passive and xPhys dimension mismatch.');

activeDomain = (P == 0);
passiveSolid = (P == 2);
activeSolid  = activeDomain & (rho >= opt.threshold);

if ~any(activeSolid(:))
    error('No active solid voxels found with threshold = %.3f.', opt.threshold);
end
fprintf('Active solid voxels = %d, passive solid voxels = %d\n', nnz(activeSolid), nnz(passiveSolid));

% Build exposed voxel surface.
fvAct = voxelSurface(activeSolid, he);
if isempty(fvAct.faces)
    error('Voxel surface extraction failed: no exposed faces were found.');
end

% Create figure.
fig = figure('Name','Active voxel topology + translucent passive (one case)', ...
    'NumberTitle','off', 'Visible','on', 'Color',opt.background, ...
    'Renderer','opengl', 'Units','pixels', 'Position',[100 80 1350 850]);
ax = axes('Parent',fig, 'Position',[0.05 0.05 0.90 0.90]);
hold(ax,'on');

% Active voxel surface with thin edges.
patch(ax, fvAct, ...
    'FaceColor', opt.activeColor, ...
    'FaceAlpha', 1.0, ...
    'EdgeColor', opt.edgeColor, ...
    'EdgeAlpha', opt.edgeAlpha, ...
    'LineWidth', opt.edgeWidth, ...
    'FaceLighting', 'flat', ...
    'AmbientStrength', opt.ambient, ...
    'DiffuseStrength', opt.diffuse, ...
    'SpecularStrength', opt.specular, ...
    'SpecularExponent', opt.specPower, ...
    'BackFaceLighting', 'lit');

% Passive region shown as a clean translucent cuboid.
if any(passiveSolid(:))
    [kz,ky,kx] = ind2sub(size(passiveSolid), find(passiveSolid));
    lower = [(min(kx)-1)*he, (min(ky)-1)*he, (min(kz)-1)*he];
    upper = [ max(kx)*he,     max(ky)*he,     max(kz)*he   ];
    verts = [lower(1) lower(2) lower(3);
             upper(1) lower(2) lower(3);
             upper(1) upper(2) lower(3);
             lower(1) upper(2) lower(3);
             lower(1) lower(2) upper(3);
             upper(1) lower(2) upper(3);
             upper(1) upper(2) upper(3);
             lower(1) upper(2) upper(3)];
    faces = [1 4 3 2; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
    patch(ax, 'Vertices', verts, 'Faces', faces, ...
        'FaceColor', opt.passiveColor, ...
        'FaceAlpha', opt.passiveAlpha, ...
        'EdgeColor', 'none', ...
        'FaceLighting', 'flat', ...
        'AmbientStrength', 0.50, ...
        'DiffuseStrength', 0.35, ...
        'SpecularStrength', 0.0, ...
        'BackFaceLighting', 'lit');
end

hold(ax,'off');
set(ax, 'Color', opt.background, 'Projection', 'orthographic', ...
    'Clipping','on', 'Visible','off');
axis(ax,'equal'); axis(ax,'vis3d'); axis(ax,'off');

L = [nx ny nz] * he;
pad = max(L) * opt.margin;
xlim(ax, [-pad, L(1)+pad]);
ylim(ax, [-pad, L(2)+pad]);
zlim(ax, [-pad, L(3)+pad]);
view(ax, opt.azimuth, opt.elevation);
camtarget(ax, L/2);
set(ax,'CameraViewAngleMode','auto','CameraPositionMode','auto');

% Camera-relative lighting so the model remains visible when rotated.
lighting(ax,'flat');
camlight(ax,'headlight');
camlight(ax,'left');

drawnow;
axis(ax,'equal');
set(ax,'CameraViewAngleMode','auto');

[folder,base] = fileparts(matFile);
outFile = fullfile(folder, [base '_ActivePassive_VoxelInteractive.fig']);
savefig(fig, outFile);
fprintf('Single-case voxel figure saved:\n%s\n', outFile);
fprintf('Figure is open. Rotate to adjust view; no PNG was generated.\n');
end

% =======================================================================
function fv = voxelSurface(mask, he)
% Extract exposed faces of a 3D voxel set as patch geometry.
% Input mask dimensions are [z, y, x].
% Output vertices use physical coordinates:
%   x along 3rd dim, y along 2nd dim, z along 1st dim.

[nz,ny,nx] = size(mask);
faces = zeros(0,4);
verts = zeros(0,3);
vertCount = 0;

% Local face patterns for a unit cube with vertices:
% 1:(x0,y0,z0), 2:(x1,y0,z0), 3:(x1,y1,z0), 4:(x0,y1,z0),
% 5:(x0,y0,z1), 6:(x1,y0,z1), 7:(x1,y1,z1), 8:(x0,y1,z1)
facePattern = {
    [1 4 8 5], ... % x- left
    [2 3 7 6], ... % x+ right
    [1 2 6 5], ... % y- front
    [4 3 7 8], ... % y+ back
    [1 2 3 4], ... % z- bottom
    [5 6 7 8]  ... % z+ top
    };

for iz = 1:nz
    z0 = (iz-1)*he; z1 = iz*he;
    for iy = 1:ny
        y0 = (iy-1)*he; y1 = iy*he;
        for ix = 1:nx
            if ~mask(iz,iy,ix)
                continue;
            end

            neigh = false(1,6);
            neigh(1) = ix>1  && mask(iz,iy,ix-1); % x-
            neigh(2) = ix<nx && mask(iz,iy,ix+1); % x+
            neigh(3) = iy>1  && mask(iz,iy-1,ix); % y-
            neigh(4) = iy<ny && mask(iz,iy+1,ix); % y+
            neigh(5) = iz>1  && mask(iz-1,iy,ix); % z-
            neigh(6) = iz<nz && mask(iz+1,iy,ix); % z+

            if all(neigh)
                continue; % interior voxel, nothing exposed
            end

            x0 = (ix-1)*he; x1 = ix*he;
            Vloc = [x0 y0 z0;
                    x1 y0 z0;
                    x1 y1 z0;
                    x0 y1 z0;
                    x0 y0 z1;
                    x1 y0 z1;
                    x1 y1 z1;
                    x0 y1 z1];
            verts = [verts; Vloc]; %#ok<AGROW>
            base = vertCount;
            vertCount = vertCount + 8;

            for f = 1:6
                if ~neigh(f)
                    faces = [faces; base + facePattern{f}]; %#ok<AGROW>
                end
            end
        end
    end
end

fv.faces = faces;
fv.vertices = verts;
end
