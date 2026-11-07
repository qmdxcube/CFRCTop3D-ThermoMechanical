function plotHistory(objectiveHistory, volumefrac, contentfracHis, objectiveLabel)
%PLOTHISTORY Plot objective and resource-constraint histories.
%   objectiveLabel is optional. Examples:
%     'Compliance'
%     'Elastic Strain Energy'

if nargin < 4 || isempty(objectiveLabel)
    objectiveLabel = 'Objective';
end

objectiveHistory = objectiveHistory(:);
volumefrac = volumefrac(:);
contentfracHis = contentfracHis(:);

n = min([numel(objectiveHistory),numel(volumefrac),numel(contentfracHis)]);
objectiveHistory = objectiveHistory(1:n);
volumefrac = volumefrac(1:n);
contentfracHis = contentfracHis(1:n);
iter = (1:n)';

figure('Color','w','Name','Iteration History');

% Left axis: optimization objective.
yyaxis left;
h1 = plot(iter,objectiveHistory,'r-','LineWidth',2);
ylabel(objectiveLabel,'FontName','Times New Roman');
ax = gca;
ax.YAxis(1).Color = 'r';

% Right axis: material and fiber volume fractions.
yyaxis right;
h2 = plot(iter,volumefrac,'b-','LineWidth',2);
hold on;
h3 = plot(iter,contentfracHis,'b--','LineWidth',2);
ylabel('Volume Fractions','FontName','Times New Roman');
ylim([0 1]);
ax.YAxis(2).Color = 'b';

% Common formatting.
grid on;
box on;
xlabel('Optimization Step','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',10);
legend([h1,h2,h3], ...
    {objectiveLabel,'Volume Fraction (Material)','Volume Fraction (Fiber)'}, ...
    'Location','northeast','Interpreter','none');
end
