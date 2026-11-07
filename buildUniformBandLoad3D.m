function [Fm, loadnid, loaddofs, loadWeights, loadInfo] = buildUniformBandLoad3D(nelx,nely,nelz,nDof,FmSum,bandHalfWidth,centerI,kIndex)
% Build a uniformly distributed z-direction mechanical load over a narrow band.
%
% By default, the original load line is
%     k = 1, i = round(nelx/2)+1, j = 1:nely+1.
% Supply centerI and kIndex when the original line is on another node column
% or surface. The load band is formed by extending that line in the x-node
% direction, and every node in the band carries the same nodal force.

if nargin < 6 || isempty(bandHalfWidth)
    bandHalfWidth = 1;
end
if nargin < 7 || isempty(centerI)
    centerI = round(nelx/2) + 1;
end
if nargin < 8 || isempty(kIndex)
    kIndex = 1;
end

bandHalfWidth = max(0,round(bandHalfWidth));
centerI = min(max(1,round(centerI)),nelx+1);
kIndex = min(max(1,round(kIndex)),nelz+1);

iMin = max(1,centerI-bandHalfWidth);
iMax = min(nelx+1,centerI+bandHalfWidth);

[kLoad,jLoad,iLoad] = ndgrid(kIndex,1:nely+1,iMin:iMax);
loadnid = sub2ind([nelz+1,nely+1,nelx+1],kLoad(:),jLoad(:),iLoad(:));
loaddofs = 3*loadnid(:);

nLoad = numel(loadnid);
loadWeights = ones(nLoad,1)/nLoad;
Fm = sparse(loaddofs,1,FmSum*loadWeights,nDof,1);

loadInfo = struct();
loadInfo.type = 'uniform-band';
loadInfo.directionDof = 'z';
loadInfo.totalLoad = FmSum;
loadInfo.bandHalfWidth = bandHalfWidth;
loadInfo.centerI = centerI;
loadInfo.midI = centerI;
loadInfo.kIndex = kIndex;
loadInfo.iRange = [iMin,iMax];
loadInfo.loadnid = loadnid;
loadInfo.loaddofs = loaddofs;
loadInfo.loadWeights = loadWeights;
end
