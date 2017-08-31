% Data generator with a control over linear separability betweeen classes
% EXAMPLES to get started
%
% For interactive demo, set options.showLevel to 1 or 2.
% -------------------------------------------------------------------------
% Version: 1.0; 2017-08-30
% Author: Nejc Ilc (nejc.ilc_at_gmail.com)

clear all; close all; clc;

%% Example 1
N = [250, 250]; % number of points in a class
K = 2; % number of classes
shapes = {'JCURVES', 'H'}; % Select shapes
% Additional options
options = [];
options.linNonSepDesiredFlag = 0; % 0: clusters must be lin. separable, 1: force lin. non. sep clusters; if ratio is [], at least one pair of clusters should be lin. non-separable
options.linNonSepDesiredAmount = 0; % parameter L: pairs of classes that are not LS (write as a negative integer)
options.minBoundDist = 0.5; % restricted zone around boundary points
options.distribution = 'UNIFORM'; % UNIFORM, GAUSS or MIXED (randomly choose between first and second)
options.randomShapesFlag = 0; % can shapes be picked on random from the list shapes?
options.showLevel = 0; % 0 - no show, 1 - only end of iter, 2 - show all
options.linprogImpl = 'MATLAB'; % 'MATLAB' (requires Optimization toolbox) or 'GLPK'

fprintf(1,'-- Example 1\nClasses: %d, min. dist.: %f, lin. non-sep. pairs: %d, distribution: %s\n', ... 
K, options.minBoundDist, options.linNonSepDesiredAmount, options.distribution);

% Create data
[data,labels,exitflag] = createDataset(K,N,shapes,options);

% Plot dataset if success
if exitflag ~= 1
    fprintf(1, 'Sorry, failed to converge.\n');
else
    pplk_scatterPlot(data,labels);
    title('Example 1');
    axis('equal');
end



%% Example 2

N = [83, 83, 83, 83, 84, 84]; % number of points in a class
K = 6; % number of classes
shapes = getShapesList('ALL'); % Select shapes

% Additional options
options = [];
options.linNonSepDesiredFlag = 0; % 0: clusters must be lin. separable, 1: force lin. non. sep clusters; if ratio is [], at least one pair of clusters should be lin. non-separable
options.linNonSepDesiredAmount = 0; % parameter L: pairs of classes that are not LS (write as a negative integer)
options.minBoundDist = 0.3; % restricted zone around boundary points
options.distribution = 'UNIFORM'; % UNIFORM, GAUSS or MIXED (randomly choose between first and second)
options.randomShapesFlag = 1; % can shapes be picked on random from the list shapes?
options.showLevel = 0; % 0 - no show, 1 - only end of iter, 2 - show all
options.linprogImpl = 'MATLAB'; % 'MATLAB' (requires Optimization toolbox) or 'GLPK'

fprintf(1,'-- Example 2\nClasses: %d, min. dist.: %f, lin. non-sep. pairs: %d, distribution: %s\n', ... 
K, options.minBoundDist, options.linNonSepDesiredAmount, options.distribution);

% Create data
[data,labels,exitflag] = createDataset(K,N,shapes,options);

% Plot dataset if success
if exitflag ~= 1
    fprintf(1, 'Sorry, failed to converge.\n');
else
    pplk_scatterPlot(data,labels);
    title('Example 2');
    axis('equal');
end


%% Example 3

N = [125, 125, 125, 125]; % number of points in a class
K = 4; % number of classes
shapes = {'J','CS','JCURVES','WAVE'}; % Select shapes

% Additional options
options = [];
options.linNonSepDesiredFlag = 1; % 0: clusters must be lin. separable, 1: force lin. non. sep clusters; if ratio is [], at least one pair of clusters should be lin. non-separable
options.linNonSepDesiredAmount = -1; % parameter L: pairs of classes that are not LS (write as a negative integer)
options.minBoundDist = 0.1; % restricted zone around boundary points
options.distribution = 'GAUSS'; % UNIFORM, GAUSS or MIXED (randomly choose between first and second)
options.randomShapesFlag = 0; % can shapes be picked on random from the list shapes?
options.showLevel = 0; % 0 - no show, 1 - only end of iter, 2 - show all
options.linprogImpl = 'MATLAB'; % 'MATLAB' (requires Optimization toolbox) or 'GLPK'

fprintf(1,'-- Example 3\nClasses: %d, min. dist.: %f, lin. non-sep. pairs: %d, distribution: %s\n', ... 
K, options.minBoundDist, options.linNonSepDesiredAmount, options.distribution);
% Create data
[data,labels,exitflag] = createDataset(K,N,shapes,options);

% Plot dataset if success
if exitflag ~= 1
    fprintf(1, 'Sorry, failed to converge.\n');
else
    pplk_scatterPlot(data,labels);
    title('Example 3');
    axis('equal');
end


%% Example 4

N = [62, 62, 62, 62, 63, 63, 63, 63]; % number of points in a class
K = 8; % number of classes
shapes = getShapesList('ALL'); % Select shapes

% Additional options
options = [];
options.linNonSepDesiredFlag = 1; % 0: clusters must be lin. separable, 1: force lin. non. sep clusters; if ratio is [], at least one pair of clusters should be lin. non-separable
options.linNonSepDesiredAmount = -2; % parameter L: pairs of classes that are not LS (write as a negative integer)
options.minBoundDist = 0.08; % restricted zone around boundary points
options.distribution = 'MIXED'; % UNIFORM, GAUSS or MIXED (randomly choose between first and second)
options.randomShapesFlag = 1; % can shapes be picked on random from the list shapes?
options.showLevel = 0; % 0 - no show, 1 - only end of iter, 2 - show all
options.linprogImpl = 'MATLAB'; % 'MATLAB' (requires Optimization toolbox) or 'GLPK'

fprintf(1,'-- Example 4\nClasses: %d, min. dist.: %f, lin. non-sep. pairs: %d, distribution: %s\n', ... 
K, options.minBoundDist, options.linNonSepDesiredAmount, options.distribution);

% Create data
[data,labels,exitflag] = createDataset(K,N,shapes,options);

% Plot dataset if success
if exitflag ~= 1
    fprintf(1, 'Sorry, failed to converge.\n');
else
    pplk_scatterPlot(data,labels);
    title('Example 4');
    axis('equal');
end
