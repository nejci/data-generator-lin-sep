function [shp,bndPnts] = getAlphaShapesMat(data,dataInd,alpha,holeThr)
% data - matrix of data points (nans are ignored)
% dataInd - vector of split points in data: [0 N1 N2 ...], N1 is number of
% data points in first cluster, etc.
% -------------------------------------------------------------------------
% Version 1.0; 2017-08-30
% Nejc Ilc (nejc.ilc_at_gmail.com)
% -------------------------------------------------------------------------

N = size(data,1);
if ~exist('dataInd','var') || isempty(dataInd)
   dataInd = [0,N];
end
assert(dataInd(end)==N,'Inconsistency.');
numClust = numel(dataInd)-1;

if ~exist('holeThr','var') || isempty(holeThr)
    holeThr = 0.01;
end

shp = cell(1,numClust);
for i = 1:numClust
    maskData = dataInd(i)+1:dataInd(i+1);
    if all(all(isnan(data(maskData,:))))
        continue;
    end
    if ~exist('alpha','var') || isempty(alpha)
        shpTmp = alphaShape(data(maskData,:));
        alpha = criticalAlpha(shpTmp,'one-region')+0.02;
    end

    shp{i} = alphaShape(data(maskData,:), alpha, 'HoleThreshold',holeThr);
end


if nargout > 1
    bndPnts = cell(1,numClust);
    for i = 1:numClust
        if isempty(shp{i})
            continue;
        end
        bndPntsTmp = boundaryFacets(shp{i});
        bndPnts{i} = bndPntsTmp(:,1);
    end

    if numClust == 1
        bndPnts = bndPnts{1};
    end
end

if numClust == 1
    shp = shp{1};
end
