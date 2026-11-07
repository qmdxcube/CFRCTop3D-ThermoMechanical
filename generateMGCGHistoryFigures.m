%% generateMGCGHistoryFigures.m
% Generate the two paper MGCG-history figures only:
%   PNG/MechanicalMGCG_q0_q1_q2_history.png  (cantilever beam)
%   PNG/MechanicalMGCG_q0_q2_q5_history.png  (double-clamped beam)

clear;
clc;
close all force;

scriptFullPath = mfilename('fullpath');
if isempty(scriptFullPath)
    dataDir = pwd;
else
    dataDir = fileparts(scriptFullPath);
end

outputDir = fullfile(dataDir,'PNG');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end
addpath(dataDir);

meanQ = [0 1 2];
strainQ = [0 2 5];
meanFiles = cell(1,3);
strainFiles = cell(1,3);

for k = 1:3
    meanFiles{k} = findThermoFileForQLocal( ...
        dataDir,'CFRCTop3D_ThermoMech_MGPCG_MGCGHistory_*.mat',meanQ(k));
    strainFiles{k} = findThermoFileForQLocal( ...
        dataDir,'CFRCTop3D_ThermoMech_StrainEnergy_MGPCG_MGCGHistory_*.mat',strainQ(k));
end

meanPrefix = fullfile(outputDir,'MechanicalMGCG_q0_q1_q2');
plotMechanicalMGCGHistory(meanFiles,meanQ,meanPrefix);
close all force;

strainPrefix = fullfile(outputDir,'MechanicalMGCG_q0_q2_q5');
plotMechanicalMGCGHistory(strainFiles,strainQ,strainPrefix);
close all force;

expectedFiles = { ...
    [meanPrefix '_history.png'], ...
    [strainPrefix '_history.png']};

for k = 1:numel(expectedFiles)
    if ~isfile(expectedFiles{k})
        error('Expected MGCG-history figure was not created: %s',expectedFiles{k});
    end
    fprintf('Saved: %s\n',expectedFiles{k});
end

fprintf('\nBoth MGCG-history figures were generated successfully.\n');

function filePath = findThermoFileForQLocal(dataDir,pattern,qExpected)
files = dir(fullfile(dataDir,pattern));
if isempty(files)
    error('No files match: %s',fullfile(dataDir,pattern));
end
matches = {};
for i = 1:numel(files)
    candidate = fullfile(files(i).folder,files(i).name);
    qFound = readQFromMatOrName(candidate);
    if isfinite(qFound) && abs(qFound-qExpected) <= 1e-10*max(1,abs(qExpected))
        matches{end+1} = candidate; %#ok<AGROW>
    end
end
if numel(matches) ~= 1
    error('Expected exactly one file for q=%g with pattern %s; found %d.', ...
        qExpected,pattern,numel(matches));
end
filePath = matches{1};
end

function q = readQFromMatOrName(filePath)
q = NaN;
try
    S = load(filePath,'q');
    if isfield(S,'q') && isscalar(S.q) && isfinite(S.q)
        q = double(S.q);
        return;
    end
catch
end
[~,name,ext] = fileparts(filePath);
tok = regexp([name ext],'_q([-+]?\d+(?:\.\d+)?)Fm','tokens','once');
if ~isempty(tok)
    q = str2double(tok{1});
end
end
