function Cnew = clusterMove(C,startPointInd,endPoint)
% Move points in C from C(startPointInd,:) to endPoint.
% If startPointInd is 'center', the reference starting point is cluster
% mean.
% -------------------------------------------------------------------------
% Version 1.0; 2017-08-30
% Nejc Ilc (nejc.ilc_at_gmail.com)
% -------------------------------------------------------------------------

if ischar(startPointInd) && strcmpi(startPointInd,'center')
    startPoint = mean(C,1);
else
    startPoint = C(startPointInd,:);
end

% Vector of movement
R = endPoint - startPoint;
Cnew = bsxfun(@plus,C,R);
