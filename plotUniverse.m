function fig = plotUniverse( aShp, universe )
% Plot universe, where clusters of stars live.
% -------------------------------------------------------------------------
% Version 1.0; 2017-08-30
% Nejc Ilc (nejc.ilc_at_gmail.com)
% -------------------------------------------------------------------------


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
