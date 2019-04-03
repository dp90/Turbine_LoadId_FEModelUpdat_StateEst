% Select eigenfrequencies and corresponding eigenvectors 
% from singular value graph and save results to current folder.
%
% output: "identified_eigenfrequencies.mat" which contains: 
% - Phi_id  : identified modes 
% - freq_id : identified frequencies
% - ind_d   : Indices of measured DOFs in the DOF vector
% 
% Author: E.Lourens & D.J.M.Fallais
%--------------------------------------------------------------------------

% Input frequencies:

Assignment_Question1;
load('Accelerations_A.mat');
FE_modelUPDATE;

% get indices for selecting eigenvectors at the desired frequencies
indx = [...
    find(freqSeg>=freq_id(1),1,'first');...
    find(freqSeg>=freq_id(2),1,'first');...
    find(freqSeg>=freq_id(3),1,'first');...
    find(freqSeg>=freq_id(4),1,'first');...
    find(freqSeg>=freq_id(5),1,'first');...
    find(freqSeg>=freq_id(6),1,'first');...
    find(freqSeg>=freq_id(7),1,'first')];

n_d = length(indx);

% Find nodal coordinates of the measured DOFs
ind_d = zeros(1,19);
NodeNr = zeros(19,1);
NodeCo = zeros(19,4);
for Nind = 1:Nd
    ind_d(Nind) = find(S_d(Nind,:));              % Indices of measured DOFs in the DOF vector
    NodeNr(Nind,1) = floor(DOF(ind_d(Nind),1));   % Node numbers
    NodeCo(Nind,:) = Nodes(NodeNr(Nind,1),:);     % Nodal coordinate
end
direction = round(rem(S_d*DOF,1)*100);            % Direction for plotting

% Obtain mode shapes from SVD of GDD at a frequency and plot
% 
Phi_id = zeros(Nd,n_d);
% figure('Name','mode shapes','units','normalized','outerposition',[0 0 .5 .3]);

% ax = gca;
%  outerpos = ax.OuterPosition;
%  ti = ax.TightInset; 
%  left = outerpos(1) + ti(1) - 0.05 ;
%  ax_width = outerpos(3) - ti(1) - ti(3) - left;
%  ax.Position = [left 0.20 ax_width 0.6];
%  hold on
 

for ind = 1:n_d
    
    % Extract and Normalize eigenvectors
    %------------------------------------------------------------------------------------
    Phi_id(:,ind) = real(squeeze(UGlobal(:,1,indx(ind))));               % obtain mode from first line of UGlobal at specified frequency index
    Phi_id(:,ind) = Phi_id(:,ind)/norm(Phi_id(:,ind));                   % normalize identified mode with identified mode   
    Phi(ind_d,ind) = Phi(ind_d,ind)/norm(Phi(ind_d,ind));                % normalize FE_mode with FE_mode    
    Phi(ind_d,ind) = sign(Phi_id(:,ind)'*Phi(ind_d,ind))*Phi(ind_d,ind); % align sign of computed mode for plotting
    
    % Determine the contributions of the identified and computed 
    % modes at sensor locations in the measured directions only
    %-----------------------------------------------------------
%     NodeCo1 = zeros(6,4);
%     NodeCo2 = zeros(6,4);
%     for i = 1:n_d
%         if direction(i) == 1
%             NodeCo1(i,:) = NodeCo(i,:)        + 20*[0, Phi_id(i,ind), 0, 0];     % X-direction: identified
%             NodeCo2(i,:) = Nodes(NodeNr(i),:) + 20*[0, Phi(ind_d(i),ind), 0, 0]; % computed             
%         elseif direction(i) == 2
%             NodeCo1(i,:) = NodeCo(i,:)        + 20*[0, 0, Phi_id(i,ind), 0];     % Y-direction: identified
%             NodeCo2(i,:) = Nodes(NodeNr(i),:) + 20*[0, 0, Phi(ind_d(i),ind) 0];  % computed           
%         elseif direction(i) == 3
%             NodeCo1(i,:) = NodeCo(i,:)        + 20*[0, 0, 0, Phi_id(i,ind)];     % Z-direction: identified
%             NodeCo2(i,:) = Nodes(NodeNr(i),:) + 20*[0, 0, 0, Phi(ind_d(i),ind)]; % computed
%         else
%             NodeCo1(i,:) = NodeCo(i,:)        + 20*[0, 0, 0, 0];
%             NodeCo2(i,:) = Nodes(NodeNr(i),:) + 20*[0, 0, 0, 0];
%         end
%     end    
%        
%     % Plot the components of the identified and modelled modes at 
%     % the measured locations in the measured directions
%     %-----------------------------------------------------------
%                                                % new figure 
%     subplot(1,n_d+1,ind)
%     plotelem(Nodes,Elements,Types,'Numbering','off');                       % plot elements
% 
%     hold on;                                             
%     HG1 = hggroup;                                                          % create a group object as child of axis 
%     HG2 = hggroup;                                                          % only important for setting the legend          
%     plotnodes(NodeCo2,'Numbering','off','Color',[1 0 0 ],'Parent',HG1);     % Plot modes: computed   - assign to group HG1 
%     plotnodes(NodeCo1,'Numbering','off','Color',[0 0 0 ],'Parent',HG2);     % Plot modes: identified - assign to group HG2 
%     
%     hold on
%     title(['\Phi_{id' num2str(ind) '} vs. \Phi_{' num2str(ind) '}'],'FontWeight','Normal','FontSize',9);
 end
% hold off
% legend([HG1,HG2],{'computed','identified'},'Location',[0.8 0.7 0.075 0.05]);
% 
% print -djpeg Idmodeshapes.jpg -r300


%% save identified eigendata for model updating step
savefile = 'identified_eigdata.mat';
save(savefile, 'Phi_id','freq_id','ind_d');


