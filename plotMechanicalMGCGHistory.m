function summaryTable = plotMechanicalMGCGHistory(caseFiles,qValues,outPrefix)
%PLOTMECHANICALMGCGHISTORY Plot mechanical-solve MGCG work for three cases.
%   summaryTable = plotMechanicalMGCGHistory(caseFiles,qValues,outPrefix)
%
%   caseFiles: cell array with three thermo-mechanical MAT files.
%   qValues:   corresponding heat-flux values.
%
%   The mechanical-solve work is defined as the sum of the MGCG work for
%   the mechanical-state and thermal-load-displacement equations.

if nargin < 1 || isempty(caseFiles)
    caseFiles = cell(1,3);
    for k = 1:3
        caseFiles{k} = selectMatFile(sprintf( ...
            'Select thermo-mechanical MGCG-history MAT file %d of 3',k));
    end
end
if isstring(caseFiles)
    caseFiles = cellstr(caseFiles);
end
if ~iscell(caseFiles) || numel(caseFiles) ~= 3
    error('caseFiles must contain exactly three MAT-file paths.');
end

if nargin < 2 || isempty(qValues)
    qValues = zeros(1,3);
    for k = 1:3
        qValues(k) = readQValue(caseFiles{k});
    end
end
qValues = double(qValues(:)');
if numel(qValues) ~= 3 || any(~isfinite(qValues))
    error('qValues must contain exactly three finite values.');
end

if nargin < 3 || isempty(outPrefix)
    outPrefix = buildDefaultPrefix(qValues);
end

commonFontSize = 16;
lineWidth = 1.5;
markerSpacing = 20;
numberOfCases = 3;
caseLabels = cell(1,numberOfCases);
for k = 1:numberOfCases
    caseLabels{k} = ['q = ' formatQ(qValues(k))];
end

mechanicalIterationHistory = cell(numberOfCases,1);
mechanicalTimeHistory = cell(numberOfCases,1);
mechanicalColumnIndices = cell(numberOfCases,1);
numberOfOptimizationIterations = zeros(numberOfCases,1);

for iCase = 1:numberOfCases
    data = loadMGCGHistory(caseFiles{iCase});
    validateQValue(data,caseFiles{iCase},qValues(iCase));
    mechanicalColumns = identifyMechanicalColumns( ...
        data.mgcgCallLabels,size(data.mgcgIterHis,2));
    mechanicalColumnIndices{iCase} = mechanicalColumns;
    mechanicalIterationHistory{iCase} = ...
        sum(data.mgcgIterHis(:,mechanicalColumns),2);
    mechanicalTimeHistory{iCase} = ...
        sum(data.mgcgTimeHis(:,mechanicalColumns),2);
    numberOfOptimizationIterations(iCase) = size(data.mgcgIterHis,1);
end

meanMGCGIterations = zeros(numberOfCases,1);
maximumMGCGIterations = zeros(numberOfCases,1);
totalMGCGIterations = zeros(numberOfCases,1);
meanMGCGTime = zeros(numberOfCases,1);
maximumMGCGTime = zeros(numberOfCases,1);
totalMGCGTime = zeros(numberOfCases,1);
for iCase = 1:numberOfCases
    it = mechanicalIterationHistory{iCase};
    tm = mechanicalTimeHistory{iCase};
    meanMGCGIterations(iCase) = mean(it);
    maximumMGCGIterations(iCase) = max(it);
    totalMGCGIterations(iCase) = sum(it);
    meanMGCGTime(iCase) = mean(tm);
    maximumMGCGTime(iCase) = max(tm);
    totalMGCGTime(iCase) = sum(tm);
end

summaryTable = table( ...
    qValues(:),numberOfOptimizationIterations,meanMGCGIterations, ...
    maximumMGCGIterations,totalMGCGIterations,meanMGCGTime, ...
    maximumMGCGTime,totalMGCGTime, ...
    'VariableNames',{'HeatFlux_q_Wm2','OptimizationIterations', ...
    'MeanMechanicalMGCGIterations','MaximumMechanicalMGCGIterations', ...
    'TotalMechanicalMGCGIterations','MeanMechanicalMGCGTimeSeconds', ...
    'MaximumMechanicalMGCGTimeSeconds','TotalMechanicalMGCGTimeSeconds'});

fig = figure('Color','w','Units','centimeters','Position',[2 2 18 20]);
ax1Position = [0.14 0.63 0.79 0.30];
ax2Position = [0.14 0.16 0.79 0.30];
ax1 = axes('Parent',fig,'Units','normalized','Position',ax1Position);
ax2 = axes('Parent',fig,'Units','normalized','Position',ax2Position);

lineColors = lines(numberOfCases);
lineStyles = {'-','--','-.'};
markers = {'o','s','^'};
plotHandles = gobjects(numberOfCases,1);

hold(ax1,'on');
for iCase = 1:numberOfCases
    y = mechanicalIterationHistory{iCase};
    x = (1:numel(y))';
    markerIndices = unique([1:markerSpacing:numel(y),numel(y)]);
    plotHandles(iCase) = plot(ax1,x,y, ...
        'LineStyle',lineStyles{iCase},'LineWidth',lineWidth, ...
        'Color',lineColors(iCase,:),'Marker',markers{iCase}, ...
        'MarkerIndices',markerIndices,'MarkerSize',5,'MarkerFaceColor','none');
end
hold(ax1,'off');
formatAxes(ax1,'Optimization iteration','Mechanical-solve MGCG iterations',commonFontSize);
legend(ax1,plotHandles,caseLabels,'Location','best','Box','off', ...
    'FontName','Times New Roman','FontSize',commonFontSize);

hold(ax2,'on');
for iCase = 1:numberOfCases
    y = mechanicalTimeHistory{iCase};
    x = (1:numel(y))';
    markerIndices = unique([1:markerSpacing:numel(y),numel(y)]);
    plot(ax2,x,y,'LineStyle',lineStyles{iCase},'LineWidth',lineWidth, ...
        'Color',lineColors(iCase,:),'Marker',markers{iCase}, ...
        'MarkerIndices',markerIndices,'MarkerSize',5,'MarkerFaceColor','none');
end
hold(ax2,'off');
formatAxes(ax2,'Optimization iteration','Mechanical-solve MGCG time (s)',commonFontSize);

addPanelLabelAt(fig,ax1,'(a)',commonFontSize,0.515);
addPanelLabelAt(fig,ax2,'(b)',commonFontSize,0.025);

pngFile = [outPrefix '_history.png'];
figFile = [outPrefix '_history.fig'];
csvFile = [outPrefix '_summary.csv'];
matFile = [outPrefix '_plotData.mat'];
exportgraphics(fig,pngFile,'Resolution',600);
savefig(fig,figFile);
writetable(summaryTable,csvFile);
save(matFile,'qValues','caseLabels','caseFiles','mechanicalColumnIndices', ...
    'mechanicalIterationHistory','mechanicalTimeHistory','summaryTable');
end

function filePath = selectMatFile(dialogTitle)
[fileName,fileFolder] = uigetfile('*.mat',dialogTitle);
if isequal(fileName,0)
    error('No MAT file was selected for: %s.',dialogTitle);
end
filePath = fullfile(fileFolder,fileName);
end

function q = readQValue(filePath)
S = load(filePath,'q');
if isfield(S,'q') && isscalar(S.q) && isfinite(S.q)
    q = double(S.q);
    return;
end
[~,name,ext] = fileparts(filePath);
tok = regexp([name ext],'_q([-+]?\d+(?:\.\d+)?)Fm','tokens','once');
if isempty(tok)
    error('Cannot infer q from file: %s',filePath);
end
q = str2double(tok{1});
end

function data = loadMGCGHistory(filePath)
required = {'mgcgIterHis','mgcgTimeHis','mgcgCallLabels'};
data = load(filePath,required{:},'q');
for k = 1:numel(required)
    if ~isfield(data,required{k})
        error('MAT file %s does not contain %s.',filePath,required{k});
    end
end
if isempty(data.mgcgIterHis) || isempty(data.mgcgTimeHis)
    error('MGCG histories are empty in MAT file: %s.',filePath);
end
if ~isequal(size(data.mgcgIterHis),size(data.mgcgTimeHis))
    error('mgcgIterHis and mgcgTimeHis must have identical sizes in %s.',filePath);
end
end

function validateQValue(data,filePath,expectedQ)
if ~isfield(data,'q') || isempty(data.q)
    return;
end
if ~isscalar(data.q) || ~isfinite(data.q)
    warning('Invalid q variable in MAT file: %s.',filePath);
    return;
end
if abs(double(data.q)-expectedQ) > 1e-10*max(1,abs(expectedQ))
    warning('Selected MAT file has q = %.12g, expected q = %.12g:\n%s', ...
        double(data.q),expectedQ,filePath);
end
end

function cols = identifyMechanicalColumns(callLabels,nColumns)
if isstring(callLabels)
    callLabels = cellstr(callLabels);
elseif ischar(callLabels)
    callLabels = cellstr(callLabels);
end
if ~iscell(callLabels)
    error('mgcgCallLabels must be a cell array, string array, or character array.');
end
normalized = strings(numel(callLabels),1);
for k = 1:numel(callLabels)
    s = lower(string(callLabels{k}));
    normalized(k) = regexprep(s,'[^a-z0-9]','');
end
idx1 = find(normalized == "mechanicalstate" | normalized == "mechanicalequilibrium",1,'first');
idx2 = find(normalized == "thermaldisplacement" | normalized == "thermalloadedisplacement",1,'first');
if isempty(idx1) || isempty(idx2)
    if nColumns >= 3
        warning('Could not identify both mechanical-solve labels; using columns 2 and 3.');
        cols = [2 3];
    else
        error('At least three MGCG-history columns are required.');
    end
else
    cols = [idx1 idx2];
end
end

function formatAxes(ax,xLabelText,yLabelText,fontSize)
set(ax,'FontName','Times New Roman','FontSize',fontSize,'LineWidth',0.8, ...
    'Box','on','Layer','top','XGrid','on','YGrid','on','GridAlpha',0.18);
xlabel(ax,xLabelText,'FontName','Times New Roman','FontSize',fontSize);
ylabel(ax,yLabelText,'FontName','Times New Roman','FontSize',fontSize);
end

function addPanelLabelAt(fig,ax,panelLabel,fontSize,labelY)
oldUnits = ax.Units;
ax.Units = 'normalized';
pos = ax.Position;
ax.Units = oldUnits;
annotation(fig,'textbox',[pos(1),labelY,pos(3),0.028], ...
    'String',panelLabel,'Interpreter','none','LineStyle','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontName','Times New Roman','FontSize',fontSize,'FitBoxToText','off');
end

function prefix = buildDefaultPrefix(qValues)
parts = cell(1,numel(qValues));
for k = 1:numel(qValues)
    parts{k} = ['q' regexprep(formatQ(qValues(k)),'\.0$','')];
end
prefix = ['MechanicalMGCG_' strjoin(parts,'_')];
end

function s = formatQ(q)
if abs(q) < 1e-12
    s = '0';
else
    s = sprintf('%.1f',q);
end
end
