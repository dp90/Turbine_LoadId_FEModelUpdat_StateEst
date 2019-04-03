function RunOptim()
% 
% part A : 	Perform a brute force assessment of the cost function on the defined 
%			parameter space.
% part B :  Run a gradient based optimization algorithm in order to find an
%			optimal parameter set which is intended to minimize the difference
%			between the identified model and the parametrized FE-model with respect.
%			to the defined cost function.
%
% Dependencies: Before running this script you need to complete "ObjFun.m"
%
% Author: D.J.M.Fallais
% -------------------------------------------------------------------------

%% A: Brute force assessment

% define variables, axis names, and scale factors
% -----------------------------------------------
% Assignment of x1 and x2 in Assignment_Question2.m
x1 = evalin('base','x1');       % Define the parameter space for each variable     
x2 = evalin('base','x2');       

x1txt = 'E_{soil}';  % Variable name for axis labels                    
x2txt = 'E_{steel}';

% Define scale factors
sclx1 = x1(1);
sclx2 = x2(1);

% Assign scale factors to workspace - to allow acces from non-nested-function 
assignin('base','sclx1',sclx1) 
assignin('base','sclx2',sclx2)

% Scale variables - these are scaled back inside ObjFun.m!
x1 = x1./sclx1;
x2 = x2./sclx2;

% Evaluate objective function on grid 
% ------------------------------------
for i = 1:length(x1)
    for j = 1:length(x2) 
        f(i,j) = ObjFun( [x1(i),x2(j)] );     
    end
end

% Call function for plotting ObjFun on grid
%------------------------------------------
plot_brute

%% B: Gradient based optimization 

% Set up shared variables with nested output function "outfun"
history.x = [];
history.fval = [];
 
% Setup optimization: set bounds and initial guess for optimization variables
% lb: lower bounds - scaled
% ub: lower bounds - scaled
% x0: initial guesses - scaled 
% f0: initial cost function value
lb = [x1(1),x2(1)];    
ub = [x1(end),x2(end)];  
x0 = [x1(1), x2(1)];
f0 = ObjFun(x0);

% Set output function and default tolerances
options = optimset('OutputFcn',@outfun);
options.TolX   = 1e-10;
options.TolCon = 1e-8;

% fmincon - interior point method - unconstrained - bounded 
%----------------------------------------------------------
[x,fval,exitflag,output,lambda,hessian] = fmincon(@ObjFun,x0,...
    [],[],[],[],lb,ub,[],options);

% Assign output to workspace
% --------------------------
assignin('base','history',history)
assignin('base','output',output)

% Call function for plotting results
% ----------------------------------
plot_gradient

% end of main function definition
% -------------------------------------------------------------------------


%% Begin sub function definition
%-------------------------------
function stop = outfun(x,optimValues,state)
    % write intermediate results
    stop = false;  % don't stop
    switch state
        case 'init'
            % do nothing
        case 'iter'
            % during iteration:
            history.fval = [history.fval; optimValues.fval];
            history.x = [history.x; x];
        case 'done'
            % do nothing 
        otherwise
    end
end

function plot_brute()
    % make contour plot and save to visualise output of gradient based method
    h1 = figure('Name','Brute_Contour','units','normalized',...
                'outerposition',[0.05 0.05 .4 .3]);
    contour(x1*sclx1,x2*sclx2,f',15,'ShowText','off')
    colorbar;
    xlabel(x1txt)
    ylabel(x2txt)
    set(gcf,'Color',[1 1 1])
    savefig(h1,'Brute_Contour.fig')

    % also plot a surface
    h2 = figure('Name','Brute_Surf','units','normalized',...
                'outerposition',[0.05 0.05 .4 .5]);
    [X1,X2] = ndgrid(x1*sclx1,x2*sclx2);
    surfc(X1,X2,f,'FaceAlpha',0.2);
    xlabel(x1txt)
    ylabel(x2txt)
    zlabel('object functional [-]')
    set(gcf,'Color',[1 1 1])
    savefig(h2, 'Brute_Surf.fig')
end

function plot_gradient()
    % plot the results in generated contour and surface plot
    % load plots obtained in previous block
    close all
    
    openfig('Brute_Contour.fig');       
        hold on; 
        plot(x0(1)*sclx1,x0(2)*sclx2,'sb','MarkerFaceColor','b') % initial    
        plot([x0(1); history.x(:,1)].*sclx1,...
            [x0(2); history.x(:,2)].*sclx2,'sg-','MarkerFaceColor','None');    % intermediate
        plot(x(1)*sclx1,x(2)*sclx2,'sr','MarkerFaceColor','r')   % final
        hold off; 
% % Choose between saving matlab figure, or jpg
    savefig('Brute_Contour.fig')
%    print -djpeg Brutecontour.jpg -r300

    openfig('Brute_Surf.fig');
        hold on
        plot3([x0(1);history.x(:,1)].*sclx1,...
              [x0(2); history.x(:,2)].*sclx2,...
              [f0; history.fval])
        scatter3([x0(1);history.x(:,1)].*sclx1,...
                 [x0(2); history.x(:,2)].*sclx2,...
                 [f0; history.fval],...
                 'ro','filled')
        hold off
% % Choose between saving matlab figure, or jpg
    savefig('Brute_Surf.fig')
%     print -djpeg BruteSurf.jpg -r300
end

% end of subfunctions - end of function
%--------------------------------------------------------------------------
end