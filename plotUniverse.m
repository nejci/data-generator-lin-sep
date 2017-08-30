function fig = plotUniverse( aShp, universe )



fig = figure();
hold on;

if exist('universe','var') && ~isempty(universe)
    rectangle('Position',[universe(1,1),universe(1,2),universe(2,1)-universe(1,1),universe(2,2)-universe(1,2)]);
    plot(mean(universe(:,1)),mean(universe(:,2)),'xk');
end

numShapes = numel(aShp);

if ~iscell(aShp) && numShapes == 1
    aShp = {aShp};
end

for s = 1:numShapes
    plot(aShp{s});
end

axis equal;
hold off;
end

