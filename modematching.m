function [modepairs, Phi_ids, freq_ids, Phi_s, freq_s] = modematching(Phi_id,freq_id ,Phi,freq, S_d)
%
% Relate the identified modes to computed modes using a MAC based matching strategy 
% comprising the following five steps:
%
%	 1) Compute a MAC matrix for a subset of the modelled modes 
%	    and all identified modes. 
%	 2) Construct "modepairs"
%	    FOR the i'th identified mode: 
%	    a) get the  "i'th" column of the MAC matrix which relates the "i'th" 
% 	       identified mode with the modelled modes
% 	    b) select the largest MAC entry from that column. The index of the
% 	       largest entry corresponds to the mode number.
% 	    c) IF the largest MAC is above a predefined minimum MAC:
% 	       write the mode numbers of the identified and best matching modelled mode 
% 	       to the rows of "modepairs"
% 	       ELSE: no modelled mode will be paired with the i'th identified mode.
% 	 3) (OPTIONAL): In case multiple modelled modes have been matched with one 
% 	    identified mode: only keep the mode pair with the highest MAC and 
% 	    delete all other mode pairings from "modepairs".
% 	 4) Generate mode selection matrices (L_id and L) based on information in 
% 	   "modepairs".
% 	 5) Format the function output: select a consistent set of modes and frequencies
%		using the selection matrices constructed in S4).
%
% - Function settings:
% -> ThMAC (LINE:72): Minimum MAC value for matching of modes. 
%
% Function input arguments:
% -  Phi_id:  measured/identified modes - normalized
% -  freq_id: measured/identified frequencies [Hz]
% -  Phi:     All modelled modes  - normalized
% -  freq:    All modelled frequencies [Hz]
%
% Function output arguements: 
% -  modepairs: array containing the indices of paired identified and 
%               modelled modes in row wise format   
% -  Phi_ids:   sorted/matched identified modes  
% -  freq_ids:  sorted/matched identified frequencies
% -  Phi_s:     sorted/matched modelled modes
% -  freq_s:    sorted/matched modelled frequencies
%
% Authour: D.J.M.Fallais 
%--------------------------------------------------------------------------


% Indices of measured DOFs in the DOF vector
n_d = size(S_d,1);
for ind = 1:n_d
     ind_d(ind) = find(S_d(ind,:));     
end


%% Construct mode pairs "modepairs"

% 1) Compute MAC of identified and modelled modes 
%--------------------------------------------------
MAC =[];
for ind1 = 1:size(Phi_id,2)+5      % for a computed mode: take more to be sure that identified modes can be present in set
    for ind2 = 1:size(Phi_id,2)    % fill MAC row wise by comparing all identified modes with the ind2'th computed mode:
        MAC(ind1,ind2)=abs(Phi(ind_d,ind1)'*Phi_id(:,ind2))^2/(norm(Phi(ind_d,ind1))^2*norm(Phi_id(:,ind2))^2);
    end
end

assignin('base','MAC',MAC) 

% % Plot the MAC matrix 
% figure('Color',[ 1 1 1]);
% imagesc(MAC); colorbar;
% xlabel('Mode nr.: Identified')
% ylabel('Mode nr.: Modelled')


% 2) Basic mode pairing
%---------------------------------------------------------------
ThMAC = 0.8;								% Set the Treshold for allowing mode pairing

modepairs = [];                             % preallocate modepairs 
for i = 1:size(MAC,2)                       % for the number of identified modes: 
	if max(MAC(:,i))>=ThMAC                 % if the max MAC in this col is alrger than ThMAC, then
		[~, j] = max(MAC(:,i));             % find the largest entry i.e. the best matching computed mode and
		modepairs = [modepairs; [i, j]];    % write info to "modepairs": i: identified mode nr., j: corresponding model mode nr.
	end
end

% 3) OPTIONAL : Remove non uniquely matched modelled modes 
% % If needed, (un)comment from here on: 
% %---------------------------------------------------------
% MMmodes = unique(modepairs(:,2));       % mode numbers of the matched modelled modes 
% 
% for i = 1:length(MMmodes)               % for each matched modelled mode: 
% 
% 	% check whether the modelled mode has been matched more than once (non-uniquely)
%     [I] = find(modepairs(:,2)==MMmodes(i)); % find indices which correspond to currenlty investigated matched modelled mode
%    
%     % If more than one match: get the MAC values for all the non-uniquely paired set of identified-modelled modes. 
%     if length(I)>1                          % if multiple matches: 
%         MACVALS = [];                       % preallocate array for MAC values of non unique matches
%         for j = 1:length(I)                 % get all relevant MAC values
%             MACVALS(j) = MAC(modepairs(I(j),2),modepairs(I(j),1));
%         end
%         [~,II] = max(MACVALS);              % Find the index in I, which corr. the pair with maximum MAC value. 
%         JJ = setdiff(1:length(I),II);       % all mode pairs other than pair II need to be deleted from modepairs.
%         modepairs(I(JJ),:) = [];            % delete non unique assigned modes
%     end
% end


%% Construct selection matrices and format output:

% 4) Make selection matrices to select pairs of modes
%-------------------------------------------------------
L_id = zeros(size(Phi_id,2),size(modepairs,1));     % preallocate selection matrix for identified modes
L    = zeros(size(Phi,2),size(modepairs,1));        % selection matrix for computed modes
for i = 1:size(modepairs,1)
    L_id(modepairs(i,1),i) = 1;     % assemble mode selection matrix: identified modes
    L(modepairs(i,2),i)    = 1;     % assemble mode selection matrix: modelled modes
end


% 5) format output: order modes and frequencies
%-------------------------------------------------------
Phi_ids  = Phi_id*L_id;                        % selected "re-ordered" identified modes
freq_ids = freq_id(modepairs(:,1));            % and frequencies
Phi_s    = Phi*L;                              % selected "re-ordered" modelled modes
freq_s   = freq(modepairs(:,2));               % and frequencies