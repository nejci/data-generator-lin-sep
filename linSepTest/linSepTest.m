function isLinSep = linSepTest(data1,data2,impl,show)
% data1 - data matrix of first cluster [numPoints x numDimensions]
% data2 - data matrix of second cluster [numPoints x numDimensions]
% impl  - implementation: 'GLPK' or 'MATLAB', leave [] for automatic
%         if numPoints <= 100, GLPK is somewhat faster
% show  - set to 1 to visualize line between data if exists
% inLinSep - 1 if data1 is linearly separable from data2
%
% GLPKMEX or MATLAB Optimization toolbox required
% Nejc Ilc, 2015-10-02

[N1,dim] = size(data1);
[N2,dim2] = size(data2);
assert(dim==dim2,'Number of dimensions must agree between datasets.');

if ~exist('show','var') || isempty(show)
   show = 0; 
end

if ~exist('impl','var') || isempty(impl)
    if max(N1,N2) > 200
        impl = 'MATLAB';
    else
        impl = 'GLPK';
    end
end

P = [data1;data2];
labels = [ones(1,N1), 2*ones(1,N2)];

% set-up linear problem
c = zeros(dim+1,1); % objectives
A = [P.*[-ones(N1,dim);ones(N2,dim)], [ones(N1,1);-ones(N2,1)]];
b = -ones(N1+N2,1);
lb = -inf(dim+1,1);
ub =  inf(dim+1,1);

isLinSep = 0;

switch upper(impl)
    case 'GLPK'        
        % GNU Linear Programming Kit
        s = 1; % 1 = minimization, -1 = maximization
        ctype = repmat('U',N1+N2,1);
        vartype = repmat('C',dim+1,1);
        % Output all GLPK messages on workspace
        param = [];
        param.msglev=1; % 1 - only errors
        ticID = tic();
        [xmin,~,status]=glpk(c,A,b,lb,ub,ctype,vartype,s,param);
        t = toc(ticID);
        if status == 5
            isLinSep = 1;
        end
        
        
    case 'MATLAB'        
        % MATLAB Optimization Toolbox
        options = optimoptions('linprog','Algorithm','dual-simplex','Display','Off');
        ticID=tic();
        [xmin,~,status]=linprog(c,A,b,[],[],lb,ub,[],options);
        t=toc(ticID);
        if status == 1
            isLinSep = 1;
        end
    
    otherwise
        error('Wrong implementation! GLPK or MATLAB available.');
        
end

if show
    fprintf(1,'%s: time = %fs\n',impl,t);    
    fig = pplk_scatterPlot(P,labels);    
    if isLinSep
        h1 = xmin(1);
        h2 = xmin(2);
        beta = xmin(3);
    
        fprintf(1,'Solution found! h1 = %f, h2 = %f, beta = %f\n',h1,h2,beta);
        
        figure(fig);
        hold on;
        % drawing the separating line
        if(h2 ~= 0)
            k = -h1/h2;
            y0 = beta/h2;
            x0 = beta/h1;
            xMax = 3;
            line([0,x0,xMax],[y0,0,xMax*k+y0]);
        else
            % vertical
            line(-beta/h1,3);
        end
        hold off;
    else
        fprintf(1,'No solution!\n');
    end
end
