% %
% % Determination of eigenfrequencies of the structure from acceleration
% % data by means of frequency domain signal decomposition (SVD)
% %             
% % -----------------------------------------------------------------------
% % -----------------------------------------------------------------------

%% Load model and measurement data
clc 
close all hidden 
beep off
clear all
tic;     

addpath(strcat(cd,'\StaBIL-2.0'))

set(0,'DefaultFigureVisible','off');
% Run FE_model to inintialize the FE model of the turbine and calculate
% eigenfrequencies and mode shapes (using the StaBIL package). 
FE_model;
set(0,'DefaultFigureVisible','on');

% Load first set of measured accelerations
load('Accelerations_A.mat');

%% Computation of Sdd and singular values
N = length(t);          % finds number of measurement times
dt = t(2) - t(1);       % dt = time interval between to measurements
T = N*dt;               % T = total time frame measured
F = 1/dt;               % Sampling frequency [Hz]

% Below the time is segmented to get several measurements
TSeg = 10/(Omega(1)/(2*pi)) ;   % 1st estimation duration of a segment 
    if floor(T/TSeg) == T/TSeg  % check if results in integer
        CountSeg = T/TSeg;
    else                            % changes the timeframe a little
        CountSeg = ceil(T/TSeg);    % ceil instead of round, little cheat
                                    %   since else TSeg would be 3.33333
                                    %   which gives problems later on
        TSeg = T/CountSeg;
        fprintf(['variable TSeg is changed from %.2f to %d.2f'...
            '. Hence, there are %d Segments. \n'],...
            10/(Omega(1)/(2*pi)) , TSeg , CountSeg);
    end

zeropad = 7500;
NSeg = TSeg / dt + zeropad;                   % Nr of samples per segment
tSeg = (0:TSeg / dt-1)*dt;                    % Time axis
dfSeg = F / NSeg;                             % Frequency bin [Hz]
freqSeg = [0:NSeg/2-1,-NSeg/2:-1] * dfSeg;    % Frequency axis
Nd = size(d,1);

% Division of time in segments and store in available frequencies
firstInd = 1;                       
dSeg = zeros(Nd,NSeg,CountSeg);     % pre-allocation
for segment = 1:CountSeg
     lastInd = firstInd+(TSeg / dt)-1;
     dpart = [d(:,firstInd:lastInd), zeros(19,zeropad)];
     disp(size(dpart));
     dSeg(:,:,segment) = fft(dpart,[],2);                   
     firstInd = lastInd;
end
clear firstInd lastInd dpart;       % remove temp variables from workspace

SValue = zeros(Nd,NSeg/2);          % pre-allocation
UGlobal = zeros(Nd,Nd,NSeg/2);      % pre-allocation

for ind = 1:NSeg/2
    Sdd = zeros(Nd,Nd);             % pre-allocation
    for segment = 1:CountSeg
       Sdd = Sdd + dSeg(:,ind,segment)*dSeg(:,ind,segment)';
    end
    [U,S,V] = svd(Sdd);             % SVD
    % Store the singular values and vectors at each frequency
    SValue(:,ind) = diag(S);
    UGlobal(:,:,ind) = U;
end

fprintf('Rank of Spectral density matrix Sdd = %d \n', rank(Sdd))

%% Plot singular values of the spectral matrix
figure('Name','spectral matrix','units','normalized','outerposition',...
    [0 0 .5 .4]);

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1) +0.03 ;
ax_width = outerpos(3) - ti(1) - ti(3) - left;
ax.Position = [left 0.2 ax_width 0.75];
hold on

% Plot singular values 
for ind=1:Nd
    semilogy(freqSeg(1:NSeg/2),SValue(ind,:));
    hold on;
end

xlabel('Frequency [Hz]');
ylabel('Singular values');
ylim([100,10^8]);
% xlim([0,15]); % for singular value plot zoom
xlim([10,13]);% for zoom to special requency
set(gca,'yscale','log') %# to set the y-axis to logarithmic
set(gca,'XGrid','on');

%% IDENTIFIED EIGENFREQUENCIES
freq_id = [0.34,0.86,5.04,6.64,9.83,11.34,12.62];

% Plot vertical eigenfrequency lines
f_0Analyt = Omega/(2*pi);
ind =1;
    while f_0Analyt(ind) < freqSeg(NSeg/2)
        semilogy([f_0Analyt(ind) f_0Analyt(ind)],ylim,...
            'b','LineWidth',1); 
        ind = ind +1;
        hold on;
    end
clear ind

% Plot red dots around identified eigenfrequencies
freq_idFreq = zeros(length(freq_id),1);
for ind = 1:length(freq_id)
    index = find(freqSeg>=freq_id(ind),1,'first');
    freq_idFreq(ind) = SValue(1,index);
end
scatter(freq_id,freq_idFreq,80,'ro');
clear index
hold off;



% CHOOSE ONE PRINT STATEMENT
% print -djpeg singularValueplot.jpg -r300
% print -djpeg singularValueplotzoom.jpg -r300
% print -djpeg svPlotomega6.jpg -r300
% 
% 
% 
% fprintf('Finished run of file %s, runtime %.2f seconds\n', mfilename , toc)

%% plot of model with nodes
% load('Accelerations_B.mat'); % for loadidentification

% figure(); 
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.1, 0.3, 0.8]);                                                     
% plotnodes(Nodes); hold on;
% plotelem(Nodes,Elements,Types,'Numbering','off'); hold on
% selecter = S_d*DOF;
% Nodesred = zeros(19,4);
% for i =1:length(selecter)
%     for j = 1:40
%         if round(selecter(i)) == Nodes(j)
%             Nodesred(i,:,:,:) = Nodes(j,:,:,:);
%         else
%             
%         end
%     end
% end
% plotnodes(Nodesred, 'Numbering','off', 'color',[1 0 0],'Markersize' , 15);
% 
% print -djpeg plotofmnodesB.jpg -r300
%}

