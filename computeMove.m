function [yNew, dNew] = computeMove(x,y,k)
% Point x is fixed, point y is movable
% k - stiffness of spring (Hooke's law: F = k*distX2Y), k on [0,1]
% yNew - new position of y
% dNew - new distance between x and yNew
% -------------------------------------------------------------------------
% Version 1.0; 2017-08-30
% Nejc Ilc (nejc.ilc_at_gmail.com)
% -------------------------------------------------------------------------

% Hooke's law
% E is unit vector of R
%E = R ./ dist_x2y;
assert(k<=1 && k>=0,'k should be on the interval [0,1].');
% R = x-y   (vector pointing towards x)
R = x - y;
% New location of y
yNew = y + k.*R;
dNew = norm(x-yNew);
% if yNew and x collide, make no move
if dNew == 0
   yNew = y;
   dNew = norm(R);
end

end
