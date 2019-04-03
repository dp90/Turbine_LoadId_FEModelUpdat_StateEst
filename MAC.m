% Calculate and plot the global MAC matrix
eigenvectors;

MACmatrix = zeros(n_d,n_d+3);
% %  modal assurance matrix from linear algebra, eigenvectors are linear
% independend, hence they are orthogonal in their space. Hence, if you
% multiply them, you should always get zero if they are not equal, and 1 if
% they are the same
for ind1 = 1:n_d
    for ind2 = 1:n_d+4
       MACmatrix(ind1,ind2) = ((Phi_id(:,ind1))'*Phi(ind_d,ind2))^2 ...
           /(norm(Phi_id(:,ind1))^2*norm(Phi(ind_d,ind2))^2);
    end
end

% figure('Name','MAC','units','normalized','outerposition',[0.1 0.1 .5 .4]);
% imagesc(MACmatrix);
% 
% set(gca,...
%  'xticklabel',{'\Phi_1' '\Phi_2' '\Phi_3' '\Phi_4' '\Phi_5' '\Phi_6' '\Phi_7' '\Phi_8' '\Phi_9' '\Phi_{10}' '\Phi_{11}'},...
%  'yticklabel',{'\Phi_{id 1}' '\Phi_{id 2}' '\Phi_{id 3}' '\Phi_{id 4}' '\Phi_{id 5}' '\Phi_{id 6}' '\Phi_{id 7}'});
% 
% xlabel('Computed modes')
% ylabel('Identified modes')
% colorMap = [linspace(1,1,256)', linspace(1,0,256)',linspace(1,0,256)'];
% colormap(colorMap);
% 
% colorbar;
% % 
% print -djpeg MAC.jpg -r300

% % MODEMATHCING FUNCTION

[modepairs, Phi_ids, freq_ids, Phi_s, freq_s] = modematching(Phi_id,freq_id ,Phi,freqSeg, S_d);