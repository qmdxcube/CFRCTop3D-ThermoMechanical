%% batchPlotCFRCTop3DResults_Integrated_v2.m
% Unified paper-figure generator for the two benchmark problems and eight cases.
%
% Expected result MAT files in the same directory as this script:
%
% Cantilever / mean compliance
%   1) CFRCTop3D_MeanCompliance_...mat                         mechanical only
%   2) CFRCTop3D_ThermoMech_MGPCG_MGCGHistory_..._q0.0...mat thermo q = 0
%   3) CFRCTop3D_ThermoMech_MGPCG_MGCGHistory_..._q1.0...mat thermo q = 1
%   4) CFRCTop3D_ThermoMech_MGPCG_MGCGHistory_..._q2.0...mat thermo q = 2
%
% Double-clamped beam / elastic strain energy
%   5) CFRCTop3D_ElasticStrainEnergy_...mat                    mechanical only
%   6) CFRCTop3D_ThermoMech_StrainEnergy_MGPCG_MGCGHistory_..._q0.0...mat
%   7) CFRCTop3D_ThermoMech_StrainEnergy_MGPCG_MGCGHistory_..._q2.0...mat
%   8) CFRCTop3D_ThermoMech_StrainEnergy_MGPCG_MGCGHistory_..._q5.0...mat
%
% The six thermo-mechanical cases generate the field/history figures used in
% Figs. 2-4 and 8-10. The two mechanical-only cases are included in the
% four-case timing comparisons (Figs. 5 and 11). The three thermo-mechanical
% cases for each problem are also used for the MGCG-history plots (Figs. 6 and 12).

clear;
clc;
close all force;

scriptFullPath = mfilename('fullpath');
if isempty(scriptFullPath)
    dataDir = pwd;
else
    dataDir = fileparts(scriptFullPath);
end

% Save all generated PNG/FIG/CSV/MAT post-processing files in a dedicated PNG folder.
outputDir = fullfile(dataDir,'PNG');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end
pngResolution = 600;
showFigures = false;
overwriteExistingFiles = true;

addpath(dataDir);

requiredFunctions = { ...
    'plotOptimizedDesign', ...
    'plotHistory', ...
    'plotTimingComparison_fourCases', ...
    'plotMechanicalMGCGHistory', ...
    'setup_3d_style'};
for k = 1:numel(requiredFunctions)
    if exist(requiredFunctions{k},'file') ~= 2
        error('Required function not found: %s.m',requiredFunctions{k});
    end
end

oldFigureVisibility = get(groot,'DefaultFigureVisible');
visibilityCleanup = onCleanup(@() ...
    set(groot,'DefaultFigureVisible',oldFigureVisibility)); %#ok<NASGU>
if showFigures
    set(groot,'DefaultFigureVisible','on');
else
    set(groot,'DefaultFigureVisible','off');
end

%% ------------------------------------------------------------------------
% Locate the eight result files.
% -------------------------------------------------------------------------
meanMechFile = findUniqueMechanicalFile( ...
    dataDir,'CFRCTop3D_MeanCompliance_*.mat','Mean-compliance mechanical-only');

strainMechFile = findUniqueMechanicalFile( ...
    dataDir,'CFRCTop3D_ElasticStrainEnergy_*.mat','Strain-energy mechanical-only');

meanQ = [0 1 2];
strainQ = [0 2 5];

meanThermoFiles = cell(1,3);
strainThermoFiles = cell(1,3);
for k = 1:3
    meanThermoFiles{k} = findThermoFileForQ( ...
        dataDir,'CFRCTop3D_ThermoMech_MGPCG_MGCGHistory_*.mat',meanQ(k));
    strainThermoFiles{k} = findThermoFileForQ( ...
        dataDir,'CFRCTop3D_ThermoMech_StrainEnergy_MGPCG_MGCGHistory_*.mat',strainQ(k));
end

fprintf('\nEight cases located successfully.\n');
fprintf('  Cantilever mechanical : %s\n',shortName(meanMechFile));
for k = 1:3
    fprintf('  Cantilever q=%g       : %s\n',meanQ(k),shortName(meanThermoFiles{k}));
end
fprintf('  Double-clamped mech.  : %s\n',shortName(strainMechFile));
for k = 1:3
    fprintf('  Double-clamped q=%g   : %s\n',strainQ(k),shortName(strainThermoFiles{k}));
end

%% ------------------------------------------------------------------------
% Generate the six thermo-mechanical result sets.
% -------------------------------------------------------------------------
for k = 1:3
    processThermoResult(meanThermoFiles{k},'Compliance',meanQ(k), ...
        'compliance','Compliance',outputDir,pngResolution,overwriteExistingFiles);
end

for k = 1:3
    processThermoResult(strainThermoFiles{k},'StrainEnergy',strainQ(k), ...
        'strainEnergy','Elastic Strain Energy',outputDir,pngResolution,overwriteExistingFiles);
end

%% ------------------------------------------------------------------------
% Four-case timing comparisons.
% -------------------------------------------------------------------------
close all force;
plotTimingComparison_fourCases( ...
    meanMechFile,meanThermoFiles,meanQ, ...
    fullfile(outputDir,'Cantilever_timing_comparison'));
close all force;

plotTimingComparison_fourCases( ...
    strainMechFile,strainThermoFiles,strainQ, ...
    fullfile(outputDir,'Doubleclamped_timing_comparison'));
close all force;

%% ------------------------------------------------------------------------
% MGCG-history comparisons for the three thermo-mechanical cases.
% -------------------------------------------------------------------------
meanMGCGPrefix = fullfile(outputDir,'MechanicalMGCG_q0_q1_q2');
plotMechanicalMGCGHistory(meanThermoFiles,meanQ,meanMGCGPrefix);
meanMGCGPNG = [meanMGCGPrefix '_history.png'];
if ~isfile(meanMGCGPNG)
    error('Expected MGCG-history figure was not created: %s',meanMGCGPNG);
end
fprintf('Saved MGCG history: %s\n',meanMGCGPNG);
close all force;

strainMGCGPrefix = fullfile(outputDir,'MechanicalMGCG_q0_q2_q5');
plotMechanicalMGCGHistory(strainThermoFiles,strainQ,strainMGCGPrefix);
strainMGCGPNG = [strainMGCGPrefix '_history.png'];
if ~isfile(strainMGCGPNG)
    error('Expected MGCG-history figure was not created: %s',strainMGCGPNG);
end
fprintf('Saved MGCG history: %s\n',strainMGCGPNG);
close all force;

fprintf('\n============================================================\n');
fprintf('All paper figures have been generated successfully.\n');
fprintf('Output directory:\n%s\n',outputDir);
fprintf('============================================================\n');

%% ========================================================================
% Local functions
% ========================================================================

function processThermoResult(matFile,objectiveTag,qValue,objectiveField, ...
    objectiveLabel,outputDir,pngResolution,overwriteExistingFiles)

fprintf('\n------------------------------------------------------------\n');
fprintf('Processing %s, q = %.1f\n',objectiveTag,qValue);
fprintf('%s\n',shortName(matFile));
fprintf('------------------------------------------------------------\n');

required = {'xPhys','yPhys','pxPhys','pyPhys','pzPhys','U','T', ...
    'volumefrac','contentfracHis',objectiveField};
data = load(matFile);
for k = 1:numel(required)
    if ~isfield(data,required{k})
        error('MAT file %s is missing variable %s.',matFile,required{k});
    end
end

if isfield(data,'he') && isscalar(data.he) && isfinite(data.he) && data.he > 0
    he = double(data.he);
else
    error('MAT file %s does not contain a valid scalar he.',matFile);
end

close all force;
plotOptimizedDesign(he,he,he,data.xPhys,data.yPhys, ...
    data.pxPhys,data.pyPhys,data.pzPhys,data.U,data.T);
drawnow;

qString = sprintf('%.1f',qValue);
prefix = sprintf('CFRCTop3D_%s_q%s',objectiveTag,qString);

saveFigureByName('Displacement Field', ...
    fullfile(outputDir,[prefix '_Displacement.png']), ...
    pngResolution,overwriteExistingFiles);
saveFigureByName('T Field', ...
    fullfile(outputDir,[prefix '_Temperature.png']), ...
    pngResolution,overwriteExistingFiles);
saveFigureByName('Fiber Design', ...
    fullfile(outputDir,[prefix '_OrientationAndFraction.png']), ...
    pngResolution,overwriteExistingFiles);
close all force;

objectiveHistory = data.(objectiveField);
plotHistory(objectiveHistory,data.volumefrac,data.contentfracHis,objectiveLabel);
drawnow;
saveFigureByName('Iteration History', ...
    fullfile(outputDir,[prefix '_History.png']), ...
    pngResolution,overwriteExistingFiles);
close all force;

% The active/passive voxel image is optional because it is a supplementary
% diagnostic rather than one of the four main subpanels. If the helper is
% available, preserve the manuscript naming convention automatically.
if exist('plotActivePassiveVoxel','file') == 2
    if isfield(data,'Passive')
        passive = data.Passive;
    else
        passive = [];
    end
    try
        voxelFig = plotActivePassiveVoxel(data.xPhys,he,passive,'off');
        voxelFile = fullfile(outputDir,[prefix '_ActivePassive_Voxel.png']);
        saveFigureHandle(voxelFig,voxelFile,pngResolution,overwriteExistingFiles);
        close(voxelFig);
    catch ME
        warning('Active/passive voxel plot skipped for %s: %s',prefix,ME.message);
    end
end
end

function filePath = findUniqueMechanicalFile(dataDir,pattern,description)
files = dir(fullfile(dataDir,pattern));
% Exclude thermo-mechanical files defensively.
keep = ~contains({files.name},'ThermoMech','IgnoreCase',true);
files = files(keep);
if isempty(files)
    error('No MAT file found for %s using pattern %s.',description,pattern);
end
if numel(files) > 1
    [~,idx] = max([files.datenum]);
    warning('Multiple files found for %s; using the newest: %s', ...
        description,files(idx).name);
    files = files(idx);
end
filePath = fullfile(files.folder,files.name);
end

function filePath = findThermoFileForQ(dataDir,pattern,targetQ)
files = dir(fullfile(dataDir,pattern));
matched = false(size(files));
for k = 1:numel(files)
    q = extractQFromName(files(k).name);
    matched(k) = isfinite(q) && abs(q-targetQ) <= 1e-10*max(1,abs(targetQ));
end
files = files(matched);
if isempty(files)
    error('No thermo-mechanical MAT file found for q = %.12g using %s.', ...
        targetQ,pattern);
end
if numel(files) > 1
    [~,idx] = max([files.datenum]);
    warning('Multiple thermo files found for q=%g; using newest: %s', ...
        targetQ,files(idx).name);
    files = files(idx);
end
filePath = fullfile(files.folder,files.name);
end

function q = extractQFromName(fileName)
tok = regexp(fileName,'_q([-+]?\d+(?:\.\d+)?)Fm','tokens','once');
if isempty(tok)
    q = NaN;
else
    q = str2double(tok{1});
end
end

function saveFigureByName(figureName,filePath,resolution,overwrite)
figs = findall(groot,'Type','figure');
idx = find(strcmp(get(figs,'Name'),figureName),1,'first');
if isempty(idx)
    error('Figure named "%s" was not created.',figureName);
end
saveFigureHandle(figs(idx),filePath,resolution,overwrite);
end

function saveFigureHandle(fig,filePath,resolution,overwrite)
if exist(filePath,'file') && ~overwrite
    fprintf('  Existing file retained: %s\n',shortName(filePath));
    return;
end
exportgraphics(fig,filePath,'Resolution',resolution);
fprintf('  Saved: %s\n',shortName(filePath));
end

function name = shortName(filePath)
[~,n,e] = fileparts(filePath);
name = [n e];
end
