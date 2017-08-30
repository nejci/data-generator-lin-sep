function Cnew = clusterRotate(C,refPointInd,angle)
% Rotate points from C around C(refPointInd,:) by angle.
% If startPointInd is 'center', the rotation center is a cluster mean.
% If missing or empty angle, rotation is by random angle
% -------------------------------------------------------------------------
% Version 1.0; 2017-08-30
% Nejc Ilc (nejc.ilc_at_gmail.com)
% -------------------------------------------------------------------------

if ~exist('angle','var') || isempty(angle)
    angle = rand()*pi*2;
end

if ischar(refPointInd) && strcmpi(refPointInd,'center')
    refPoint = mean(C,1);
else
    refPoint = C(refPointInd,:);
end

% translate refPoint to origin
Cnew = bsxfun(@minus,C,refPoint);
% rotate
cosA = cos(angle);
sinA = sin(angle);
Cnew = [cosA -sinA; sinA cosA] * Cnew';
% translate back to original mean
Cnew = bsxfun(@plus,Cnew',refPoint);
