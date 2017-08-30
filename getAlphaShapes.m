function [shp,bndPnts] = getAlphaShapes(data,alpha,holeThr)
% -------------------------------------------------------------------------
% Version 1.0; 2017-08-30
% Nejc Ilc (nejc.ilc_at_gmail.com)
% -------------------------------------------------------------------------
numData = 1;
if iscell(data)
    numData = numel(data);
else
    data = {data};
end
if ~exist('holeThr','var') || isempty(holeThr)
    holeThr = 0.01;
end

shp = cell(1,numData);
for i = 1:numData
    if isempty(data{i})
        continue;
    end
    if ~exist('alpha','var') || isempty(alpha)
        shpTmp = alphaShape(data{i});
        alpha = criticalAlpha(shpTmp,'one-region')+0.02;
    end

    shp{i} = alphaShape(data{i}, alpha, 'HoleThreshold',holeThr);

end
if numData == 1
    shp = shp{1};
end


if nargout > 1
    bndPnts = cell(1,numData);
    for i = 1:numData
        if isempty(shp{i})
            continue;
        end
        bndPntsTmp = boundaryFacets(shp{i});
        bndPnts{i} = bndPntsTmp(:,1);
    end

    if numData == 1
        bndPnts = bndPnts{1};
    end
end
