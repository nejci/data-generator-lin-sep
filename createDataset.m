function [data,labels,exitflag,moreInfo] = createDataset(K,N,shapes,options)
% [data,labels,exitflag,moreInfo] = createDataset(K,N,shapes,options)
% Algorithm for sythetic data generation with a control over linear
% separability betweeen classes.
%
% See examples.m to get started.
%
% K - number of clusters
% N - number of points in clusters; scalar or a vector [1xK]
% shapes - cell string of shapes that can be used for clusters:
%         A, B, C(S), H, I, J, L, O(S), S, T, U(S), V, Y, Z
%         LINE(S), CROSS, BARU1(S), BARU2(S), BARO(S), JCURVE(S), HALFRING(S)
%         RECT1, RECT2, SNOWMAN, DOT(S), OVAL, WAVE(S)
% options - struct with optional parameters:
%         showLevel: [0] - no show, 1 - only end of iter, 2 - show all
%         stiffness: [0.5] how many percent of distance is reduced on each iteration
%         minBoundDist: [0.5] restricted zone around boundary points
%         minBoundDistTol: [minBoundDist * 0.1] tolerance, stopping criterion
%         strategyReferencePoints: ['nearest2center'] 'center','nearest','nearest2center'
%         angleCoarseMax: [pi] rotation on a move
%         angleFineMax: [pi/2] rotation on fine-tuning when avoiding collision
%         numTrialsMax: [3] number of fresh starts if goal is not reached
%         numItersMax: [300] max. number of iterations of move on a coarse level
%         numItersFineTuneMax: [200] max. number of iterations on a fine level (avoiding collision)
%         linNonSepDesiredFlag: [0] 0: clusters must be lin. separable, 1: force lin. non. sep clusters; if ratio is [], at least one pair of clusters should be lin. non-separable
%         linNonSepDesiredAmount: [] how many cluster pairs should be linearly non-separable? As integer (-2 means 2 pairs) or as a ratio [0,1]
%         distribution: ['UNIFORM'] UNIFORM, GAUSS or MIXED (randomly choose between first and second)
%         randomScaleAmount:  [0.5] randomly scale clusters on creation? scale factor = 1 - amount; 0 - no scale, 0.2 - scale factor on [0.8,1], 1 - scale factor on [0,1]
%         randomShapesFlag: [1] can shapes be picked on random from the list shapes?
%         randomShapesOrderFlag: [1] random order of shapes to appear in simulation?
%         linprogImpl: ['GLPK'] for small problems (N ~ 200) faster than MATLAB
% data - matrix [numPoints,2] of data points
% labels - vector of labels (1...K) for each datapoint
% exitflag - 1: success, requirements are met; -1: failure
% moreInfo - additional info
% -------------------------------------------------------------------------
% Version 1.0; 2017-08-30
% Nejc Ilc (nejc.ilc_at_gmail.com)
% -------------------------------------------------------------------------

% Options defaults:
showLevel = 0; % 0 - no show, 1 - only end of iter, 2 - show all
stiffness = 0.5; % how many percent of distance is reduced on each iteration
minBoundDist = 0.5; % restricted zone around boundary points
minBoundDistTol = minBoundDist * 0.1; % tolerance, stopping criterion
strategyReferencePoints = 'nearest2center'; % 'center','nearest','nearest2center'
angleCoarseMax = pi; % rotation on a move
angleFineMax = pi/2; % rotation on fine-tuning when avoiding collision
numTrialsMax = 3; % number of fresh starts if goal is not reached
numItersMax = 300; % max. number of iterations of move on a coarse level
numItersFineTuneMax = 200; % max. number of iterations on a fine level (avoiding collision)
linNonSepDesiredFlag = 0; % 0: clusters must be lin. separable, 1: force lin. non. sep clusters; if ratio is [], at least one pair of clusters should be lin. non-separable
linNonSepDesiredAmount = []; % how many cluster pairs should be linearly non-separable? As integer (-2 means 2 pairs) or as a ratio [0,1]
distribution = 'UNIFORM'; % UNIFORM, GAUSS or MIXED (randomly choose between first and second)
randomScaleAmount = 0.5; % randomly scale clusters on creation? scale factor = 1 - amount; 0 - no scale, 0.2 - scale factor on [0.8,1], 1 - scale factor on [0,1]
randomShapesFlag = 1; % can shapes be picked on random from the list shapes?
randomShapesOrderFlag = 1; % random order of shapes to appear in simulation?
linprogImpl = 'GLPK'; % for small problems (N ~ 200) faster than MATLAB

% User values
if(exist('options','var'))
    fldNames = fieldnames(options);
    for i=1:numel(fldNames)
        fName = fldNames{i};
        switch lower(fName)
            case 'showlevel'
                showLevel = options.(fName);
            case 'stiffness'
                stiffness = options.(fName);
            case 'minbounddist'
                minBoundDist = options.(fName);
            case 'minbounddisttol'
                minBoundDistTol = options.(fName);
            case 'strategyreferencepoints'
                strategyReferencePoints = options.(fName);
            case 'anglecoarsemax'
                angleCoarseMax = options.(fName);
            case 'anglefinemax'
                angleFineMax = options.(fName);
            case 'numtrialsmax'
                numTrialsMax = options.(fName);
            case 'numitersmax'
                numItersMax = options.(fName);
            case 'numitersfinetunemax'
                numItersFineTuneMax = options.(fName);
            case 'linnonsepdesiredflag'
                linNonSepDesiredFlag = options.(fName);
            case 'linnonsepdesiredamount'
                linNonSepDesiredAmount = options.(fName);
            case 'distribution'
                distribution = options.(fName);
            case 'randomscaleamount'
                randomScaleAmount = options.(fName);
            case 'randomshapesflag'
                randomShapesFlag = options.(fName);
            case 'randomshapesorderflag'
                randomShapesOrderFlag = options.(fName);
            case 'linprogimpl'
                linprogImpl = options.(fName);
        end
    end
end
options = [];
options.showLevel = showLevel;
options.stiffness = stiffness;
options.minBoundDist = minBoundDist;
options.minBoundDistTol = minBoundDistTol;
options.strategyReferencePoints = strategyReferencePoints;
options.angleCoarseMax = angleCoarseMax;
options.angleFineMax = angleFineMax;
options.numTrialsMax = numTrialsMax;
options.numItersMax = numItersMax;
options.numItersFineTuneMax = numItersFineTuneMax;
options.linNonSepDesiredFlag = linNonSepDesiredFlag;
options.linNonSepDesiredAmount = linNonSepDesiredAmount;
options.distribution = distribution;
options.randomScaleAmount = randomScaleAmount;
options.randomShapesFlag = randomShapesFlag;
options.randomShapesOrderFlag = randomShapesOrderFlag;
options.linprogImpl = linprogImpl;
%--------------------------------------------------------------------------


ticID = tic();

if ~exist('shapes','var') || isempty(shapes)
    shapes = getShapesList('all');
end
numShapes = numel(shapes);

% outer loop - trials (fresh restarts if exitflag == -1)
numTrials = 0;
while numTrials < numTrialsMax
    numTrials = numTrials+1;
    if showLevel > 0
        fprintf(1,'Trial %d\n',numTrials);
    end

    if randomShapesFlag == 1
        shpInd = randi(numShapes,1,K);
    else
        if numShapes < K
            rep = ceil(K/numShapes);
            shpInd = repmat(1:numShapes,1,rep);
            shpInd = shpInd(1:K);
        else
            shpInd = 1:K;
        end
    end
    if randomShapesOrderFlag == 1
        shpInd = shpInd(randperm(K));
    end

    if numel(N) == 1
        N = repmat(N,1,K);
    end

    % Universe
    universe = [0 0; 50 50];
    universeCenter = mean(universe,1);
    spawnRadius = 20;

    % Lin non sep setup
    % Prepare a matrix for evidence of lin. non. separability between
    % cluster pairs.
    linNonSepMat = false(K,K);
    numPairsAll = K*(K-1)/2;
    if ~isempty(linNonSepDesiredAmount) && linNonSepDesiredAmount < 0
        % absolute value -> compute ratio
        % bound value
        linNonSepDesiredAmount = min(abs(linNonSepDesiredAmount),numPairsAll);
        linNonSepDesiredAmount = linNonSepDesiredAmount/numPairsAll;
    end


    % dataset
    data = nan(sum(N),2);
    dataInd = [0,cumsum(N)];

    % Random scaling
    % lowerBound is defined by randomScaleAmount
    % upperBound is defined by number of samples (50 -> 250)
    assert(randomScaleAmount >= 0 && randomScaleAmount <= 1,'randomScaleAmount out of bounds [0,1].');
    lowerBound = (1-randomScaleAmount);
    upperBound = 0.5*((N(1)-50)/200) +0.5;
    scaleClusterVal =  min(1,(upperBound-lowerBound) * rand() + lowerBound);  % random between lowerBound and 1

    exitflag = 0; % 1: all conditions met, -1: num iters over

    % first cluster is placed in the center of universe
    clust = clusterCreate(N(1), shapes{shpInd(1)},'scale',scaleClusterVal,distribution);
    clust = clusterRotate(clust,'center'); % random rotation
    clust = clusterMove(clust,'center',universeCenter);
    data(dataInd(1)+1:dataInd(2),:) = clust;

    numItersHist = ones(1,K);

    % second and further clusters
    % - spawn the cluster on the ring around universe center
    % - apply Hooke's law - force that moves second cluster towards the
    % first one and is proportional to distance between them
    for k = 2:K

        upperBound = 0.5*((N(k)-50)/200) +0.5;
        scaleClusterVal =  min(1,(upperBound-lowerBound) * rand() + lowerBound);
        clust = clusterCreate(N(k), shapes{shpInd(k)},'scale',scaleClusterVal,distribution);
        % Position the cluster randomly on the ring around the center of universe
        randAngle = rand()*pi*2;
        spawnPoint = [cos(randAngle), sin(randAngle)] .* spawnRadius + universeCenter;

        clust = clusterMove(clust,'center',spawnPoint);
        clust = clusterRotate(clust,'center'); % random rotation
        data(dataInd(k)+1:dataInd(k+1),:) = clust;

        % get alphaShapes with border points (ind) from clusters
        [aShp,bndPoints] = getAlphaShapesMat(data,dataInd);

        % Simulation
        simulationRunning = 1;
        numIters = 0;
        while simulationRunning
            numIters = numIters +1;

            if showLevel == 2
                if exist('fig','var')
                    close(fig);
                end
                fig = plotUniverse(aShp,universe);
                input(['Start of iteration ',num2str(numIters)]);
            end

            %------------------------------------------------------------------
            % Choose reference points
            % a) choose one boundary point from fixed clusters and one from moving cluster
            % b) cluster center
            % c) nearest points
            fixInd = randi(k-1);
            movInd = k;

            numBndPoints = cellfun(@numel,bndPoints);

            switch strategyReferencePoints
                case 'random'
                    % random boundary points from fix and mov
                    numBndPointsFix = numBndPoints(fixInd);
                    numBndPointsMov = numBndPoints(movInd);
                    refPointFixInd = bndPoints{fixInd}(randi(numBndPointsFix));
                    refPointMovInd = bndPoints{movInd}(randi(numBndPointsMov));
                    refPointFix = aShp{fixInd}.Points(refPointFixInd,:);
                    refPointMov = aShp{movInd}.Points(refPointMovInd,:);

                case 'center'
                    refPointFix = mean(aShp{fixInd}.Points,1); % mean of boundary points?
                    refPointMov = mean(aShp{movInd}.Points,1); % mean of boundary points?

                case 'nearest'
                    % nearest boundary points between fix and mov cluster
                    [nI,distMov2Fix] = nearestNeighbor(aShp{fixInd}, aShp{movInd}.Points(bndPoints{movInd},:));
                    [~,refPointMovInd] = min(distMov2Fix);
                    refPointFixInd = nI(refPointMovInd);
                    refPointMovInd = bndPoints{movInd}(refPointMovInd);
                    refPointFix = aShp{fixInd}.Points(refPointFixInd,:);
                    refPointMov = aShp{movInd}.Points(refPointMovInd,:);

                case 'nearest2center'
                    % nearest boundary point of mov to center of fix boundary
                    % Why center of boundary? Faster. Also, with this strategy
                    % we want po penetrate to inside of fix cluster and mean of
                    % all points is not always the center of shape.
                    [~,distMov2Fix] = nearestNeighbor(aShp{fixInd}, aShp{movInd}.Points(bndPoints{movInd},:));
                    [~,refPointMovInd] = min(distMov2Fix);
                    refPointMovInd = bndPoints{movInd}(refPointMovInd);
                    refPointFix = mean(aShp{fixInd}.Points(bndPoints{fixInd},:),1);
                    refPointMov = aShp{movInd}.Points(refPointMovInd,:);

                otherwise
                    error('Wrong strategyReferencePoints');
            end

            %------------------------------------------------------------------
            % Try to move cluster (moving) towards fixed one, randomly rotate
            % Compute direction and amount of movement of moving cluster

            CmovOld = aShp{movInd}.Points; % restore this state if new move fails

            [movNew,distNew] = computeMove(refPointFix,refPointMov,stiffness);
            angleCoarse = angleCoarseMax*(rand()*2-1) * distNew/spawnRadius;
            if strcmpi(strategyReferencePoints,'center')
                CmovNew = clusterMove(aShp{movInd}.Points,'center',movNew);
                CmovNew = clusterRotate(CmovNew,'center',angleCoarse);
            else
                CmovNew = clusterMove(aShp{movInd}.Points,refPointMovInd,movNew);
                CmovNew = clusterRotate(CmovNew,refPointMovInd,angleCoarse);
            end

            if showLevel == 2
                close(fig);
                fig = plotUniverse(aShp);
                figure(fig);
                hold on;
                plot(refPointFix(1),refPointFix(2),'gx');
                plot(refPointMov(1),refPointMov(2),'bx');
                input('Choosen points ...');

                close(fig);
                aShp{movInd}.Points = CmovNew; % update alpha shape points
                fig = plotUniverse(aShp);
                input(['Try to make this move: distNew = ',num2str(distNew),', angleCoarse = ',num2str(angleCoarse*180/pi),'�']);
            end


            %------------------------------------------------------------------
            % Check for contraints about distance between boundary points
            numItersFineTune = 0;
            isFineTuned = 0;
            while ~isFineTuned

                numItersFineTune = numItersFineTune+1;

                % gather boundary points of fixed clusters
                boundPointsCmp = zeros(sum(numBndPoints(1:movInd-1)),2);
                bndPntsInd = [0 cumsum(numBndPoints)];
                for s = 1:movInd-1
                    boundPointsCmp(bndPntsInd(s)+1:bndPntsInd(s+1),:) = aShp{s}.Points(bndPoints{s},:);
                end
                %isOverlap2 = inShape(aShp{movInd},boundPointsCmp(:,1),boundPointsCmp(:,2));

                % boundary points of mov cluster
                boundPointsMov = CmovNew(bndPoints{movInd},:);

                % Compute minimal distance between boundary points mov:others
                dMov2Oth = sqrt(sqdistance2(boundPointsMov,boundPointsCmp));
                [minValTmp,minInd] = min(dMov2Oth,[],1);
                [minVal,minIndFix] = min(minValTmp);
                minIndMov = minInd(minIndFix);
                CminIndMov = bndPoints{movInd}(minIndMov);

                % If minimal distance is lower than threshold, move cluster
                % away:
                % a) in the direction of the nearest boundary point
                % b) in the direction of the center of mov cluster
                diffBndDist = minVal - minBoundDist;
                if diffBndDist < 0
                    % move point backwards
                    pntFix = boundPointsCmp(minIndFix,:);
                    pntMov = boundPointsMov(minIndMov,:);
                    vecMove = (pntMov - pntFix)/minVal; % unit vector
                    pntMovNew = pntMov + vecMove*abs(diffBndDist);
                    CmovNew = clusterMove(CmovNew,CminIndMov,pntMovNew);
                    % compute angle of fine rotation: larger is the penetration
                    % inside restricted zone, larger is the rotation and vise
                    % versa
                    %angleFine = angleFineMax*(rand()*2-1) * abs(diffBndDist)/minBoundDist;
                    angleFine = angleFineMax*...
                        (rand()*2-1) * ... % ((-1)^round(rand()))
                        (1-numItersFineTune/numItersFineTuneMax); % weight by number of iterations
                    %abs(diffBndDist)/minBoundDist * ... % weight by closeness to boundary


                    CmovNew = clusterRotate(CmovNew,'center',angleFine);

                    if showLevel == 2
                        close(fig);
                        aShp{movInd}.Points = CmovNew;
                        fig = plotUniverse(aShp);
                        input(['Boundary points too close, moved back: distNew = ',...
                            num2str(norm(pntFix-pntMovNew)),...
                            ', angleFine = ', num2str(angleFine*180/pi),'�']);
                    end

                else
                    isFineTuned = 1;
                end

                % watchdog for fine-tunning
                if numItersFineTune > numItersFineTuneMax
                    if showLevel > 0
                        warning(['numItersFineTune exceeded max value (',num2str(numItersFineTuneMax),'); not fine-tuned.']);
                    end
                    break;
                end
            end

            % update of alpha shape
            aShp{movInd}.Points = CmovNew;
            % check overlap between mov cluster and fixed data
            isOverlap = any(inShape(aShp{movInd},data(1:dataInd(k),:)));

            % check minBoundDist for every point in mov cluster
            isOverlapDeep = 0;
            if ~isOverlap && isFineTuned
                Dmat = sqrt(sqdistance2(aShp{movInd}.Points,data(1:dataInd(k),:)));
                diffBndDist = min(Dmat(:)) - minBoundDist;

                if diffBndDist < -minBoundDistTol
                    isOverlapDeep = 1;
                end
            end

            if isOverlap || isOverlapDeep || ~isFineTuned
                % reset move if clusters overlap or if fine-tunning is not
                % resolved
                aShp{movInd}.Points = CmovOld;
                data(dataInd(k)+1:dataInd(k+1),:) = CmovOld;
                if showLevel > 0
                    fprintf(1,'\tOverlap or not fine-tuned, reseting move.\n');
                end
            else
                % Update data with mov cluster new position
                data(dataInd(k)+1:dataInd(k+1),:) = CmovNew;
            end


            %------------------------------------------------------------------
            % Check stopping criteria

            % watchdog for number of iterations
            stopCrit_numIters = 0;
            if numIters > numItersMax
                if showLevel > 0
                    warning(['numIters exceeded max value (',num2str(numItersMax),'), exiting simulation for k=',num2str(k)]);
                end
                stopCrit_numIters = 1;
            end

            % if minimal distance between boundary points is exactly as desired
            stopCrit_boundDist = 0;
            if isFineTuned && (abs(diffBndDist) <= minBoundDistTol)
                %fprintf(1,'Stopping criterion met: minBoundDistTol reached.\n');
                stopCrit_boundDist = 1;
            end

            % Linear seperability test between clusters boundaries
            % For all the pairs of mov->fix clusters, run linSepTest
            % Compute only for mov cluster.
            oldPath=chdir('linSepTest');
            for fixI=1:k-1
                linNonSepMat(movInd,fixI) = linSepTest(...
                    aShp{movInd}.Points(bndPoints{movInd},:),...
                    aShp{fixI}.Points(bndPoints{fixI},:),linprogImpl) == 0;
            end
            chdir(oldPath);
            linNonSepNum = sum(linNonSepMat(:));
            linNonSepRatio = linNonSepNum/numPairsAll; % global ratio
            isAnyLinNonSep = linNonSepNum > 0;

            stopCrit_linNonSep = 0;
            if linNonSepDesiredFlag == 0 || linNonSepDesiredAmount == 0
                if ~isAnyLinNonSep
                    stopCrit_linNonSep = 1;
                end
            else
                if isAnyLinNonSep
                    % check if current mov cluster already made its best
                    if any(linNonSepMat(movInd,:))
                        stopCrit_linNonSep = 1;
                    end

                    if isempty(linNonSepDesiredAmount)
                        % is there at least one nonSep pair?
                        stopCrit_linNonSep = 1;
                    else
                        % do we met desired ratio of nonSep pairs?
                        if linNonSepRatio == linNonSepDesiredAmount
                            stopCrit_linNonSep = 1;
                        end
                        % desired amount of non-linearity is exceeded
                        if linNonSepRatio > linNonSepDesiredAmount
                            stopCrit_linNonSep = 0;
                            % reset move
                            aShp{movInd}.Points = CmovOld;
                            data(dataInd(k)+1:dataInd(k+1),:) = CmovOld;
                            % change strategy
                            strategyReferencePoints = 'nearest';
                        end
                    end
                end
            end

            % aggregate stopping criteria
            if stopCrit_numIters || (stopCrit_linNonSep && stopCrit_boundDist)
                simulationRunning = 0;
                % exitflag
                if stopCrit_numIters || ...
                    (linNonSepDesiredFlag && ~isempty(linNonSepDesiredAmount) && ...
                    linNonSepRatio < linNonSepDesiredAmount)
                    exitflag = -1;
                else
                    exitflag = 1;
                end
            end

            % show
            if showLevel > 0
                if exist('fig','var')
                    close(fig);
                end
                fig = plotUniverse(aShp);
                fprintf(1,'End of iteration %d, fine-tune iters: %d\n',numIters,numItersFineTune);
                fprintf(1,'simulationRunning: %d, linNonSepNum: %d, linNonSepRatio: %f\n',simulationRunning,linNonSepNum,linNonSepRatio);
                input('Press enter ...');
            end
        end
        numItersHist(k) = numIters;
    end

    if exitflag == 1
        break;
    end
end

% scale on [0,1] proportionaly
data = normalizePropor(data);

% create labels vector
labels = zeros(size(data,1),1);
for k=1:K
    labels(dataInd(k)+1:dataInd(k+1)) = k;
end

moreInfo = [];
moreInfo.K = K;
moreInfo.N = N;
moreInfo.options = options;
moreInfo.linNonSepRatio = linNonSepRatio;
moreInfo.minBoundDistReached = minVal;
moreInfo.numItersHist = numItersHist;
moreInfo.numTrials = numTrials;
moreInfo.shapesList = shapes;
moreInfo.shapesUsed = shapes(shpInd);
moreInfo.timeElapsed = toc(ticID);



% output message
if showLevel > 0
    fprintf(1,'Dataset created with ');
    if exitflag == 1
        fprintf(1,'[SUCCESS]\n');
    else
        fprintf(1,'[FAIL] exitflag: %d\n',exitflag);
    end

    fprintf('\tReached linNonSepRatio: %f (flag: %d, desired: %f)\n',linNonSepRatio,linNonSepDesiredFlag,linNonSepDesiredAmount);
    fprintf('\tReached minBoundDist: %f (desired: %f, tol:%f)\n',minVal,minBoundDist, minBoundDistTol);
    fprintf('\tIterations: %d (allowed: %d)\n',numIters,numItersMax);

    pplk_scatterPlot(data,labels);
    axis equal;
end

% if showLevel > 0
%
%     % labelsClsKM = pplk_runClusterer('KM',data,K,1);
%     % labelsClsSPECLS = pplk_runClusterer('SPECLS',data,K,1);
%     % opt=[];
%     % opt.axisStyle = 'equal';
%     % pplk_scatterPlot(data,[labelsClsKM,labelsClsSPECLS],[],opt);
% end
