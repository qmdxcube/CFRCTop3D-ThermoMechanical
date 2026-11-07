function [stageSummary,caseSummary] = plotTimingComparison_fourCases( ...
    mechFile,thermoFiles,qValues,outPrefix)
%PLOTTIMINGCOMPARISON_FOURCASES Compare four optimization cases.
%   [stageSummary,caseSummary] = plotTimingComparison_fourCases( ...
%       mechFile,thermoFiles,qValues,outPrefix)
%
%   thermoFiles must contain exactly three thermo-mechanical MAT files.
%   qValues must contain the corresponding three heat-flux values.
%
%   The output figure contains:
%     (a) average computational time per optimization iteration;
%     (b) accumulated computational time over the optimization process.

if nargin < 1 || isempty(mechFile)
    mechFile = selectMatFile('Select the mechanical-only result MAT file');
end
if nargin < 2 || isempty(thermoFiles)
    thermoFiles = cell(1,3);
    for k = 1:3
        thermoFiles{k} = selectMatFile(sprintf( ...
            'Select thermo-mechanical MAT file %d of 3',k));
    end
end
if isstring(thermoFiles)
    thermoFiles = cellstr(thermoFiles);
end
if ~iscell(thermoFiles) || numel(thermoFiles) ~= 3
    error('thermoFiles must contain exactly three MAT-file paths.');
end

if nargin < 3 || isempty(qValues)
    qValues = zeros(1,3);
    for k = 1:3
        qValues(k) = readQValue(thermoFiles{k});
    end
end
qValues = double(qValues(:)');
if numel(qValues) ~= 3 || any(~isfinite(qValues))
    error('qValues must contain exactly three finite values.');
end

if nargin < 4 || isempty(outPrefix)
    outPrefix = 'CFRCTop3D_timing_comparison';
end

commonFontSize = 12;
stageNames = { ...
    'Filtering/projection', ...
    'Material evaluation', ...
    'Matrix assembly', ...
    'Thermal solve', ...
    'Mechanical solve', ...
    'Sensitivity analysis', ...
    'MMA update', ...
    'Orientation update'};

mechTime = loadTimingHistory(mechFile,7,'Mechanical-only');
mechAligned = [mechTime(:,1:3),zeros(size(mechTime,1),1),mechTime(:,4:7)];

thermoAligned = cell(1,3);
for k = 1:3
    thermoAligned{k} = loadTimingHistory(thermoFiles{k},8, ...
        sprintf('Thermo-mechanical q = %.12g',qValues(k)));
end

averageTimes = [ ...
    mean(mechAligned,1); ...
    mean(thermoAligned{1},1); ...
    mean(thermoAligned{2},1); ...
    mean(thermoAligned{3},1)];

accumulatedTimes = [ ...
    sum(mechAligned,1); ...
    sum(thermoAligned{1},1); ...
    sum(thermoAligned{2},1); ...
    sum(thermoAligned{3},1)];

iterationCounts = [ ...
    size(mechAligned,1); ...
    size(thermoAligned{1},1); ...
    size(thermoAligned{2},1); ...
    size(thermoAligned{3},1)];

averageTotalTimes = sum(averageTimes,2);
accumulatedTotalTimes = sum(accumulatedTimes,2);
averageOverheadPercent = 100*(averageTotalTimes-averageTotalTimes(1))/ ...
    max(averageTotalTimes(1),eps);
accumulatedOverheadPercent = 100*(accumulatedTotalTimes-accumulatedTotalTimes(1))/ ...
    max(accumulatedTotalTimes(1),eps);

caseNames = [ ...
    "Mechanical only"; ...
    "Thermo-mechanical, q = " + string(formatQ(qValues(1))); ...
    "Thermo-mechanical, q = " + string(formatQ(qValues(2))); ...
    "Thermo-mechanical, q = " + string(formatQ(qValues(3)))];

stageSummary = table( ...
    string(stageNames(:)), ...
    averageTimes(1,:)',averageTimes(2,:)',averageTimes(3,:)',averageTimes(4,:)', ...
    accumulatedTimes(1,:)',accumulatedTimes(2,:)',accumulatedTimes(3,:)',accumulatedTimes(4,:)', ...
    'VariableNames',{ ...
    'Stage','MechanicalAverageSecondsPerIteration', ...
    'ThermoCase1AverageSecondsPerIteration','ThermoCase2AverageSecondsPerIteration', ...
    'ThermoCase3AverageSecondsPerIteration','MechanicalAccumulatedSeconds', ...
    'ThermoCase1AccumulatedSeconds','ThermoCase2AccumulatedSeconds', ...
    'ThermoCase3AccumulatedSeconds'});

caseSummary = table( ...
    caseNames,iterationCounts,averageTotalTimes,averageOverheadPercent, ...
    accumulatedTotalTimes,accumulatedOverheadPercent, ...
    'VariableNames',{'Case','NumberOfIterations','AverageSecondsPerIteration', ...
    'AverageOverheadRelativeToMechanicalPercent','AccumulatedIterationSeconds', ...
    'AccumulatedOverheadRelativeToMechanicalPercent'});

writetable(stageSummary,[outPrefix '_stage_timing_summary.csv']);
writetable(caseSummary,[outPrefix '_case_timing_summary.csv']);

fig = figure('Color','w','Units','centimeters','Position',[2 0.5 27 24]);
ax1Position = [0.10 0.61 0.67 0.32];
ax2Position = [0.10 0.14 0.67 0.32];
stackColors = lines(numel(stageNames));

caseLabels = cell(1,4);
caseLabels{1} = {'Mechanical only'};
for k = 1:3
    caseLabels{k+1} = {'Thermo-mechanical',['q = ' formatQ(qValues(k))]};
end

ax1 = axes('Parent',fig,'Units','normalized','Position',ax1Position);
b1 = bar(ax1,averageTimes,'stacked','BarWidth',0.76);
formatBarSeries(b1,stackColors);
formatAxes(ax1,4,'Average time per iteration (s)',commonFontSize);
addTotalLabels(ax1,averageTotalTimes,commonFontSize);

ax2 = axes('Parent',fig,'Units','normalized','Position',ax2Position);
b2 = bar(ax2,accumulatedTimes,'stacked','BarWidth',0.76);
formatBarSeries(b2,stackColors);
formatAxes(ax2,4,'Accumulated iteration time (s)',commonFontSize);
addTotalLabels(ax2,accumulatedTotalTimes,commonFontSize);

lgd = legend(ax2,b2,stageNames,'Location','none','Box','off', ...
    'FontName','Times New Roman','FontSize',commonFontSize);
lgd.Title.String = 'Computational stage';
lgd.Title.FontSize = commonFontSize;
lgd.Title.FontName = 'Times New Roman';
drawnow;
ax1.Position = ax1Position;
ax2.Position = ax2Position;
lgd.Units = 'normalized';
p = lgd.Position;
p(1) = 0.80;
p(2) = ax2Position(2) + 0.5*(ax2Position(4)-p(4));
lgd.Position = p;

addCaseLabelsBelow(fig,ax1,caseLabels,commonFontSize);
addCaseLabelsBelow(fig,ax2,caseLabels,commonFontSize);
addPanelLabelBelow(fig,ax1,'(a)',commonFontSize);
addPanelLabelBelow(fig,ax2,'(b)',commonFontSize);

pngFile = [outPrefix '_four_cases_stacked_timing.png'];
figFile = [outPrefix '_four_cases_stacked_timing.fig'];
exportgraphics(fig,pngFile,'Resolution',600);
savefig(fig,figFile);
end

function filePath = selectMatFile(dialogTitle)
[name,path] = uigetfile('*.mat',dialogTitle);
if isequal(name,0)
    error('No MAT file was selected for: %s.',dialogTitle);
end
filePath = fullfile(path,name);
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

function timeHis = loadTimingHistory(filePath,expectedColumns,caseDescription)
S = load(filePath,'timeHis');
if ~isfield(S,'timeHis') || isempty(S.timeHis)
    error('%s MAT file does not contain a nonempty timeHis.',caseDescription);
end
timeHis = double(S.timeHis);
lastUsed = find(any(timeHis > 0,2),1,'last');
if isempty(lastUsed)
    error('%s timeHis contains no positive timing values.',caseDescription);
end
timeHis = timeHis(1:lastUsed,:);
if size(timeHis,2) ~= expectedColumns
    error('%s timeHis must contain %d columns, but contains %d.', ...
        caseDescription,expectedColumns,size(timeHis,2));
end
if any(~isfinite(timeHis),'all') || any(timeHis < 0,'all')
    error('%s timeHis contains invalid values.',caseDescription);
end
end

function formatBarSeries(b,colorMap)
for k = 1:numel(b)
    b(k).FaceColor = colorMap(k,:);
    b(k).EdgeColor = 'none';
end
end

function formatAxes(ax,nCases,yLabelText,fontSize)
ax.XTick = 1:nCases;
ax.XTickLabel = repmat({''},1,nCases);
ax.XTickLabelRotation = 0;
ax.FontName = 'Times New Roman';
ax.FontSize = fontSize;
ax.LineWidth = 0.8;
ax.Box = 'on';
ax.Layer = 'top';
ax.YGrid = 'on';
ax.XGrid = 'off';
ax.GridAlpha = 0.18;
ax.XLim = [0.45,nCases+0.55];
ylabel(ax,yLabelText,'FontName','Times New Roman','FontSize',fontSize);
end

function addCaseLabelsBelow(fig,ax,caseLabels,fontSize)
oldUnits = ax.Units;
ax.Units = 'normalized';
pos = ax.Position;
ax.Units = oldUnits;
xLimits = ax.XLim;
spacing = pos(3)/(xLimits(2)-xLimits(1));
boxWidth = 0.92*spacing;
boxHeight = 0.052;
boxY = pos(2)-0.060;
for k = 1:numel(caseLabels)
    relativeX = (k-xLimits(1))/(xLimits(2)-xLimits(1));
    xCenter = pos(1)+relativeX*pos(3);
    annotation(fig,'textbox',[xCenter-boxWidth/2,boxY,boxWidth,boxHeight], ...
        'String',caseLabels{k},'Interpreter','none','LineStyle','none', ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontName','Times New Roman','FontSize',fontSize,'FitBoxToText','off');
end
end

function addPanelLabelBelow(fig,ax,label,fontSize)
oldUnits = ax.Units;
ax.Units = 'normalized';
pos = ax.Position;
ax.Units = oldUnits;
y = pos(2)-0.115;
annotation(fig,'textbox',[pos(1),y,pos(3),0.028], ...
    'String',label,'Interpreter','none','LineStyle','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontName','Times New Roman','FontSize',fontSize,'FitBoxToText','off');
end

function addTotalLabels(ax,totals,fontSize)
yMax = max(totals);
for k = 1:numel(totals)
    text(ax,k,totals(k)+0.02*max(yMax,eps),sprintf('%.2f s',totals(k)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontName','Times New Roman','FontSize',fontSize);
end
ylim(ax,[0,1.12*max(yMax,eps)]);
end

function s = formatQ(q)
if abs(q) < 1e-12
    s = '0';
else
    s = sprintf('%.1f',q);
end
end
