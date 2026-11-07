clc;
clear;
close all;

%% RunAll
% Sequentially runs the eight benchmark cases used for the first two
% numerical examples:
%
% Problem 1: Mean compliance
%   Case 1 - Mechanical-only loading
%   Case 2 - Thermo-mechanical loading, q = 0 W/m^2
%   Case 3 - Thermo-mechanical loading, q = 1 W/m^2
%   Case 4 - Thermo-mechanical loading, q = 2 W/m^2
%
% Problem 2: Elastic strain energy
%   Case 5 - Mechanical-only loading
%   Case 6 - Thermo-mechanical loading, q = 0 W/m^2
%   Case 7 - Thermo-mechanical loading, q = 2 W/m^2
%   Case 8 - Thermo-mechanical loading, q = 5 W/m^2
%
% The q = 0 thermo-mechanical cases are retained as consistency checks
% against the corresponding mechanical-only formulations.

%% Common parameters
lx = 1.0;
nelx = 80;
nely = 40;
nelz = 40;
he = lx/nelx;

volfrac = 0.30;
contentfrac = 0.09;

rxmin = 1.5;
rymin = 1.5;
rpmin = 1.5;

cgtol = 1e-5;
cgmax = 1000;
maxiter = 200;

FmSum = -1.0e6;
loadBandHalfWidth = 3.0;

totalTimer = tic;

%% Problem 1: Mean compliance
fprintf('\n============================================================\n');
fprintf('Case 1/8: Mean Compliance - mechanical only\n');
fprintf('============================================================\n\n');
CFRCTop3D_MeanCompliance( ...
    nelx,nely,nelz,he,volfrac,contentfrac, ...
    rxmin,rymin,rpmin,cgtol,cgmax,maxiter,FmSum);
close all;

qMean = [0.0, 1.0, 2.0];
for iq = 1:numel(qMean)
    q = qMean(iq);
    fprintf('\n============================================================\n');
    fprintf('Case %d/8: Mean Compliance - thermo-mechanical, q = %.1f W/m^2\n', ...
        iq+1,q);
    fprintf('============================================================\n\n');
    CFRCTop3D_ThermoMech_MeanCompliance( ...
        nelx,nely,nelz,he,volfrac,contentfrac, ...
        rxmin,rymin,rpmin,cgtol,cgmax,maxiter,q,FmSum);
    close all;
end

%% Problem 2: Elastic strain energy
fprintf('\n============================================================\n');
fprintf('Case 5/8: Elastic Strain Energy - mechanical only\n');
fprintf('============================================================\n\n');
CFRCTop3D_ElasticStrainEnergy( ...
    nelx,nely,nelz,he,volfrac,contentfrac, ...
    rxmin,rymin,rpmin,cgtol,cgmax,maxiter,FmSum,loadBandHalfWidth);
close all;

qStrainEnergy = [0.0, 2.0, 5.0];
for iq = 1:numel(qStrainEnergy)
    q = qStrainEnergy(iq);
    fprintf('\n============================================================\n');
    fprintf('Case %d/8: Elastic Strain Energy - thermo-mechanical, q = %.1f W/m^2\n', ...
        iq+5,q);
    fprintf('============================================================\n\n');
    CFRCTop3D_ThermoMech_StrainEnergy( ...
        nelx,nely,nelz,he,volfrac,contentfrac, ...
        rxmin,rymin,rpmin,cgtol,cgmax,maxiter,q,FmSum,loadBandHalfWidth);
    close all;
end

fprintf('\n============================================================\n');
fprintf('All eight benchmark cases completed.\n');
fprintf('Total wall-clock time: %.3f s\n',toc(totalTimer));
fprintf('============================================================\n');
