function list = getShapesList(mode)
% list = getShapesList(mode)
% mode: ['all'], 'low', 'medium', 'high' - ability to embed another cluster
% -------------------------------------------------------------------------
% Version 1.0; 2017-08-30
% Nejc Ilc (nejc.ilc_at_gmail.com)
% -------------------------------------------------------------------------

if ~exist('mode','var') || isempty(mode)
    mode = 'all';
end

switch lower(mode)
    case 'all'
        list = {'A', 'B', 'C','CS', 'H', 'I', 'J', 'L', 'O','OS', 'S', 'T', ...
            'U', 'US', 'V', 'Y', 'Z', 'LINE', 'LINES', 'CROSS', 'BARU1','BARU1S', ...
            'BARU2','BARU2S','BARO','BAROS','JCURVE','JCURVES','HALFRING','HALFRINGS',...
            'RECT1','RECT2','SNOWMAN','DOT','DOTS','OVAL','WAVE','WAVES'};

    % shapes selected for PHD thesis
    case 'phd_demo'
        list = {'A', 'CS', 'H', 'OS', 'LINES', 'CROSS','SNOWMAN','DOTS','WAVE'};

    % groups based on probability of embedding another cluster
    case 'low'
        list = {'I', 'LINE', 'LINES', 'DOT', 'DOTS', 'RECT1', 'RECT2','SNOWMAN'};

    case 'medium'
        list = {'B','C','H','J','O','S','T','Y','Z',...
            'BARO','BARU1','BARU2','HALFRING','JCURVE','OVAL','WAVE',};

    case 'high'
        list = {'A','CS','L','OS','U','US','V',...
            'CROSS','BARU1S','BARU2S','BAROS','JCURVES',...
            'HALFRINGS','WAVES'};
    case 'slim'
        list = {'CS','OS','US','LINES','BARU1S','BARU2S','BAROS','JCURVES',...
            'HALFRINGS','DOTS','WAVES'};
    case 'sphere'
        list = {'RECT1','RECT2','DOT','DOTS','O','OS','OVAL','BARO','BAROS'};
    case 'longitudinal'
        list = {'RECT2','I','J','L','T','Y','LINE','LINES','BARU1','BARU1S','JCURVE','JCURVES',...
        'HALFRING','HALFRINGS','SNOWMAN','WAVE','WAVES'};
end
