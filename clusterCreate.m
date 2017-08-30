function data = clusterCreate(N,shape,trans,transVal,distribution,show)
% data = clusterCreate(N,shape,trans,transVal,distribution,show)
% Creates set of points, i.e. a cluster.
%
% N: number of points
% shape: shape of cluster
%         A, B, C(S), H, I, J, L, O(S), S, T, U(S), V, Y, Z  
%         LINE(S), CROSS, BARU1(S), BARU2(S), BARO(S), JCURVE(S), HALFRING(S)
%         RECT1, RECT2, SNOWMAN, DOT(S), OVAL, WAVE(S)
% trans: transformation type
%        'none'
%        'mirror' - mirror cluster around X and Y-axis
%        'rotate' - rotate cluster
%        'scale'  - scale down cluster size
% transVal: value of transformation
% distribution: 'UNIFORM', 'GAUSS' or 'MIXED' to randomly choose among the first two
% show: 1 - verbatim on, 0 - off
%--------------------------------------------------------------------------
% Nejc Ilc, 2015

if ~exist('shape','var') || isempty(shape)
    shapesList = getShapesList('all');
    shape = shapesList{randi(numel(shapesList))};
end

if ~exist('trans','var') || isempty(trans)
    trans = 'none';
end
if ~exist('distribution','var') || isempty(distribution)
    distribution = 'uniform';
end
if ~exist('show','var') || isempty(show)
    show = 0;
end

data = zeros(N,2);

blockSize = N*2;

mirrorY = 0;
rescale = 1; % rescale proportionaly on [0,1]

if strcmpi(distribution,'MIXED')
   % choose UNIFORM or GAUSS
   c = randi(2);
   if c == 1
       distribution = 'UNIFORM';
   else
       distribution = 'GAUSS';
   end
end



i = 0;
numIters = 0;
while i < N
    numIters = numIters+1;
    
    % distribution
    switch upper(distribution)
        case 'UNIFORM'
            x = rand(blockSize,1);
            y = rand(blockSize,1);
        case 'GAUSS'
            x = 0.4.*randn(blockSize*2,1)+rand();
            y = 0.4.*randn(blockSize*2,1)+rand();
            maskX=(~(0 <= x & x <= 1));
            maskY=(~(0 <= y & y <= 1));
            out = maskX | maskY;                        
            % erase outliers
            x(out)=[];
            y(out)=[];
            
        otherwise
            error('Wrong distribution name.');
    end
    switch upper(shape)
        case 'A'
            mask = (y >= 1 - 2*x) & (y >= 2*x - 1) & ...
                ( ((y <= 1.3 - 2*x) | (y <= 2*x - 0.7)) | (y > 0.5 & y <= 0.7) );
            mirrorY = 1;            
        case 'B'
            up = sqrt( (x-0.5).*(x-0.5) + (y-0.3).*(y-0.3) );
            down = sqrt( (x-0.5).*(x-0.5) + (y-0.7).*(y-0.7) );
            mask = (x <= 0.2) | ...
                (x < 0.5 & y <= 0.15) | ...
                (x < 0.5 & y >= 0.85) | ...
                (x < 0.5 & y > 0.4 & y < 0.6) | ...
                (x >= 0.5 & up <= 0.3 & up >= 0.15 ) | ...
                (x >= 0.5 & down <= 0.3 & down >= 0.15);
            mirrorY = 1;
        case 'C'
            dist = sqrt( (x-0.5).*(x-0.5) + (y-0.5).*(y-0.5) );
            mask = (x <= 0.5) & (dist <= 0.5) & (dist >= 0.3) | ...
                ((x > 0.5) & (x < 0.9) & (y < 0.3 | y >= 0.7) & ...
                (dist <= 0.5) & (dist >= 0.3) );
            mirrorY = 1;
        case 'CS'
            dist = sqrt( (x-0.5).*(x-0.5) + (y-0.5).*(y-0.5) );
            mask = (x <= 0.5) & (dist <= 0.5) & (dist >= 0.4) | ...
                ((x > 0.5) & (x < 0.9) & (y < 0.3 | y >= 0.7) & ...
                (dist <= 0.5) & (dist >= 0.4) );
            mirrorY = 1;
        case 'H'
            mask = (x <= 0.2) | (x >= 0.7 & x <= 0.9) | ...
                (y >= 0.4 & y < 0.6 & x <= 0.9);
            mirrorY = 1;
        case 'I'
            mask = (x >= 0.4 & x <= 0.6) | ...
                (y < 0.2 & x > 0.3 & x < 0.7) | ...
                (y > 0.8 & x > 0.3 & x < 0.7);
            mirrorY = 1;
        case 'J'
            dist = sqrt( (x-0.3).*(x-0.3) + (y-0.7).*(y-0.7) );
            
            mask = (x > 0.4 & x < 0.6 & y <= 0.7) | ...
                (y < 0.2 & x >= 0.3 & x <= 0.7) | ...
                (dist <= 0.3  & dist > 0.1 & y > 0.7);
            mirrorY = 1;
        case 'L'
            mask = (x < 0.2) | (x <0.8 & y > 0.8);
            mirrorY = 1;
        case 'O'
            dist = sqrt( (x-0.5).*(x-0.5) + (y-0.5).*(y-0.5) );
            mask = (dist < 0.5) & (dist >= 0.3);
            mirrorY = 1;  
        case 'OS'
            dist = sqrt( (x-0.5).*(x-0.5) + (y-0.5).*(y-0.5) );
            mask = (dist < 0.5) & (dist >= 0.4);
            mirrorY = 1;
        case 'S'
            up = sqrt( (x-0.5).*(x-0.5) + (y-0.3).*(y-0.3) );
            down = sqrt( (x-0.5).*(x-0.5) + (y-0.7).*(y-0.7) );
            
            mask = (up <= 0.3 & up > 0.1 & ...
                (x <= 0.5 | (x > 0.5 & y < 0.3))) | ...
                (down <= 0.3 & down > 0.1 & ...
                (x >= 0.5 | (x < 0.5 & y > 0.7)));
            mirrorY = 1;        
        case 'T'
            mask = (x > 0.4 & x < 0.6) | ...
                (y < 0.2 & x > 0.1 & x < 0.9);
            mirrorY = 1;            
        case 'U'
            dist = sqrt( (x-0.5).*(x-0.5) + (y-0.5).*(y-0.5) );
            mask = (y < 0.7  & (x < 0.2 | x > 0.8)) | ...
                (y>= 0.7 & (dist <= 0.5 & dist > 0.3));
            mirrorY = 1;   
        case 'US'
            dist = sqrt( (x-0.5).*(x-0.5) + (y-0.5).*(y-0.5) );
            mask = (y < 0.7  & (x < 0.1 | x > 0.9)) | ...
                (y>= 0.7 & (dist <= 0.5 & dist > 0.4));
            mirrorY = 1; 
        case 'V'
            mask = (y <= 2*x & y <= -2*x + 2)  &...
                (y >2*x-0.3 | y > -2*x + 1.7);
            mirrorY = 1;            
        case 'Y'
            mask = (y < 0.4 & y >= 0.8-x & y <= 1-x) |...
                (y < 0.4 & y >= x-0.2 & y <=x) |...
                (y >=0.4 & x > 0.4 & x < 0.6);
            mirrorY = 1;            
        case 'Z'
            mask = (x <= 0.8) & ((y < 0.2) | (y > 0.8) | (y > 0.8-x & y <= 1-x));
            mirrorY = 1;            
        case 'CROSS' % +
            mask = (x >= 0.4 & x <= 0.6) |...
                (y >= 0.4 & y <= 0.6);            
        case 'BARU1' % U-bar with short wings
            mask = (x <= 0.2 & y <=0.4) | y <= 0.2 | (x >= 0.8 & y <= 0.4);
        case 'BARU1S' % slim
            mask = (x <= 0.1 & y <=0.4) | y <= 0.1 | (x >= 0.9 & y <= 0.4);            
        case 'BARU2' % U-bar with long wings
            mask = x <= 0.2 | y <= 0.2 | x >= 0.8;
        case 'BARU2S' % slim
            mask = x <= 0.1 | y <= 0.1 | x >= 0.9;            
        case 'BARO' % O-bar
            mask = (x < 0.2) | (x > 0.8) | ...
		     (y < 0.2) | (y > 0.8);         
		case 'BAROS' % O-bar slim
            mask = (x < 0.1) | (x > 0.9) | ...
		     (y < 0.1) | (y > 0.9) ;         
        case 'LINE' % Straight line, 0.2 thick
            mask = y <= 0.6 & y >= 0.4;
        case 'LINES' % Straight line, 0.1 thick 
            mask = y <= 0.55 & y >= 0.45;
        case 'RECT1' % square           
            mask = true(size(x,1),1);        
        case 'RECT2' % Rectangle 1 x 0.6       
            mask = y <= 0.8 & y >= 0.2;            
        case 'SNOWMAN' % Two circles
            dist1 = sqrt((x-0.5).*(x-0.5) + (y-0.2).*(y-0.2));
            dist2 = sqrt((x-0.5).*(x-0.5) + (y-0.6).*(y-0.6));		
            mask = dist1 < 0.2 | dist2 < 0.3 | ...
                (x >= 0.5 & x < 0.6 & y <=0.9) ;            
        case 'DOT' % Full circle
            dist = sqrt((x-0.5).*(x-0.5) + (y-0.5).*(y-0.5));
            mask = dist <= 0.5;        
        case 'DOTS' % Full circle
            dist = sqrt((x-0.5).*(x-0.5) + (y-0.5).*(y-0.5));
            mask = dist <= 0.5;
            s = 0.25;
            x = x.*s+0.5-0.5*s;
            y = y.*s+0.5-0.5*s;
            rescale = 0;            
        case 'OVAL'           
            dist = sqrt((x-0.5).*(x-0.5) + (y-0.5).*(y-0.5));
            mask = dist <= 0.5;
            y = y*2;
        case 'JCURVE' % J
            dist = sqrt((x-0.5).*(x-0.5) + (y-0.5).*(y-0.5)); 
            mask = ((x < 0.5) & (dist <= 0.5) & (dist >= 0.3)) |...
			 ((x >= 0.5) & (y <= 0.3-0.1));
        case 'JCURVES'
            dist = sqrt((x-0.5).*(x-0.5) + (y-0.5).*(y-0.5)); 
            mask = ((x < 0.5) & (dist <= 0.5) & (dist >= 0.4)) |...
			 ((x >= 0.5) & (y <= 0.1));         
        case 'HALFRING' % Halfring, 0.2 thick
            dist = sqrt((x-0.5).*(x-0.5) + (y-0.5).*(y-0.5));            
            mask = x < 0.6 & dist <= 0.5 & dist > 0.3;
        case 'HALFRINGS' % Halfring, 0.1 thick
            dist = sqrt((x-0.5).*(x-0.5) + (y-0.5).*(y-0.5));            
            mask = x < 0.6 & dist <= 0.5 & dist > 0.4; 
        case 'WAVE'
            sine = 0.1.*sin(3*pi.*x)+0.5;
            mask = (y > sine-0.05) & (y < sine+0.05);            
            rescale = 0;
        case 'WAVES'
            sine = 0.1.*sin(3*pi.*x)+0.5;
            mask = (y > sine-0.025) & (y < sine+0.025);            
            rescale = 0;    
        otherwise
            error('Wrong shape.');
            
            
            
    end
    newPoints = [x(mask), y(mask)];
    
    % remove duplicates
    newPoints = unique(newPoints,'rows');
    % 
    count = size(newPoints,1);
    newPointsNum = min(N - i, count);
    % shuffle newPoints
    newPoints = newPoints(randperm(count),:);
    data(i+1:i+newPointsNum,:) = newPoints(1:newPointsNum,:);
    i = i+newPointsNum;
end



assert(N==size(data,1),'Not enough points!');

switch lower(trans)
    case 'mirror'
        if ~exist('transVal','var') || isempty(transVal)
            transVal = randi(3)-1;
        end        
        if transVal == 0            
        elseif transVal == 1
            data(:,1) = 1-data(:,1);
        elseif transVal == 2
            data(:,2) = 1-data(:,2);
        else
            error('Wrong transVal for mirror.');
        end
        
    case 'rotate'
        if ~exist('transVal','var') || isempty(transVal)
            transVal = rand()*2*pi;
        end
        
        % convert to radians
        %transValRad = -transVal * pi /180;
        
        % compute mean
        center = mean(data,1);
        % translate mean to origin
        data = bsxfun(@minus,data,center);
        % rotate
        cosA = cos(transVal);
        sinA = -sin(transVal);
        data = [cosA -sinA; sinA cosA] * data';
        % translate back to original mean
        data = bsxfun(@plus,data',center);
        
    case 'scale'
        if ~exist('transVal','var') || isempty(transVal)
            transVal = rand();
        end
        if transVal > 1 || transVal < 0
           error('Scaling factor only on [0,1].'); 
        end
        % compute mean
        center = mean(data,1);
        % translate mean to origin
        data = bsxfun(@minus,data,center);
        % scale (shrink)
        data = data.*transVal;
        % translate back to original mean
        data = bsxfun(@plus,data,center);
        % turn off rescaling on [0,1]
        rescale = 0;        
end

if mirrorY
    % mirror over y-axis once more (implementation trick)
    data(:,2) = 1-data(:,2);
end

if rescale
    % rescale on (0,1) interval
    data = normalizePropor(data);
end

if show
    fprintf(1,'Iterations: %d\n',numIters);
    pplk_scatterPlot(data);
    axis('equal');
end
