function [f] = ObjFun(x)
%
% Evaluate the objective/cost function at the parameter value(s) contained in x.
% 
% Required completions:
%
% 	Formulate a cost function which will be used in "RunOptim.m"; This 
% 	cost function has to expressess the similarity between a set of selected
% 	identified frequencies/modes and the corresponding modelled frequencies/modes. 
%
% 	You are free to form a cost function which is different to the provided template.
% 	For the construction of the cost function you could use the following ingredients:
%
% 	-> The identified and computed modes and frequencies without any pre processing: 
% 	- freq_id(:,i) : identified i^th frequency [1 x ni]     ni: nr. of identified modes and frequencies
% 	- Phi_id(i)    : identified i^th mode shape [nd x ni]   nd: nr. of sensors
% 	- freq(i)      : computed i^th frequency [ndof]
% 	- Phi(:,i)     : computed i^th mode shape [ndof x ndof]
%
% 	-> A set of identified and modelled modes and frequencies which are matched using
% 	a MAC based matching procedure, which is provided in modematching.m: 
% 	- freq_ids(i)  : matched identified i^th frequency [1 x ni]
% 	- Phi_ids(:,i) : matched identified i^th mode [nd x ni]	
% 	- freq_s(i)    : matched computed i^th frequency [ndof]
% 	- Phi_s(:,i)   : matched computed i^th mode [ndof x ndof]
%
% Function input: 
% - x: vector containing variables for optimisation
%
% Function output: 
% - f: cost function value
%
% Author: D.J.M.Fallais
% -------------------------------------------------------------------------

% Get scaling variables from workspace
sclx1 = evalin('base','sclx1');
sclx2 = evalin('base','sclx2');

% Load saved results from identification (as in practical 4) 
% - Phi_id(:,i) : identified i^th mode  
% - freq_id(i)  : identified i^th frequency
% - ind_d       : index
load identified_eigdata.mat

% Calculate Phi and Omega for set of design variables using FE-model
% call "[M,K,DOF,Omega,Phi] = FE_fun(X);" to evaluate model at current parameter set: 
% output, only use Phi, Omega:
[~,~,~,Omega,Phi] = FE_fun([x(1)*sclx1, x(2)*sclx2]);  
freq = Omega/(2*pi);    

assignin('base','OmegaFEM',Omega);  % write to workspace to be able to recover after optimization
assignin('base','PhiFEM',Phi);      % write to workspace to be able to recover after optimization

% Pair calculated modes with identified modes using a MAC based matching:
S_d = evalin('base','S_d');																			% load the measurement selection matrix from the WS 
[modepairs, Phi_ids, freq_ids, Phi_s, freq_s] = modematching(Phi_id, freq_id, Phi, freq, S_d);			% Match/pair identified modes with computed modes

%% Construct and compute the cost function value:
nmatch = size(modepairs,1);     																		% Number of matched modes
maxnomod = 5;																							% Define a maxium number of modes to be used

for i = 1:min(maxnomod,nmatch)     % for the number of matched modes or less.
    
	% Normalize computed modes and correct for potential sign switching 
    Phi_s2(:,i) = (S_d*Phi_s(:,i))/norm(S_d*Phi_s(:,i)); 	 % Normalize full computed mode (Phi_id has been scaled in eigenvectors). 
    Phi_s2(:,i) = sign(Phi_ids(:,i)'*Phi_s2(:,i))*Phi_s2(:,i);  % If opposite signs, switch sign of computed mode: based on inner product between identified and computed

	% define a cost function contributions related to e.g. the frequency and/or the mode shapes:
    T1(i) = (norm(Phi_s2(:,i) - Phi_ids(:,i))^2) / (norm(Phi_ids(:,i))^2);
    T2(i) = (freq_s(i) - freq_ids(i))^2 / (freq_ids(i)^2);   

end

% Sum Cost function contributions:
% In order to keep the cost function smooth, it could be important to account
% for a possibly varying number of matched modes.
f = sum(T1) + sum(T2);









