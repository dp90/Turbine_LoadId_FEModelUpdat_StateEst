% % Load identification with least square, with and without Tikhonov
% % regularization. 
% %
% % -----------------------------------------------------------------------
% % -----------------------------------------------------------------------

clc 
close all hidden 
beep off
clear variables

set(0,'DefaultFigureVisible','off');
load('Accelerations_B.mat');
FE_modelUPDATE; % use the updated finite element model
load('Forces.mat');
set(0,'DefaultFigureVisible','on');

%% FRF 
zeta = 0.02;  % Constant modal damping ratio assumed
Cstar = diag(2.*Omega.*zeta);
C = M*Phi*Cstar*Phi'*M;

c15=selectdof(DOF,15.01);   % Fwave
c35=selectdof(DOF,35.01);   % Fwind
S_p = [c15;c35].';          % [np x nDOF]

p = S_p*[Fwave;Fwind];
n_p = 2;
pw = fft(p,[],2);

% % Sampling parameters
 N = size(t,2); 
 N2 = N/2-1;
 dt = t(2)-t(1);
 F = 1/dt;
 df = F/N;
 freqHz =[0:N2,-N/2:-1]*df;
 freqRad = freqHz*2*pi;
 H = zeros(105,105,N/2);              % tranfer function to 105 x 105 DOF
 A = zeros(5,2,N/2);                 % tranfer function to 19 x 2 DOF
 %%
 dw = fft(d,[],2);                 
 
 for ind = 2:N/2
    H(:,:,ind) = (-freqRad(ind).^2)*inv(-freqRad(ind)^2*M+1i*freqRad(ind)*C+K);  
    A(:,:,ind) = S_d * squeeze(H(:,:,ind)) * S_p;
 end

% % END construction FRF or TRANSFER F.
%% Plot of FRF
%
%  figure('Name','FRF','units','normalized','outerposition',...
%         [0 0 .5 .4]);
% % subplot(2,1,1)
% title('Amplitude')
% plot(freqRad(1:N/2),abs(squeeze(A(1,1,:))))
% % subplot(2,1,2)
% % plot(freqRad(1:N/2),angle(squeeze(A(1,1,:))))
% % ylim(0 -pi)
% % yticks([-0.5*pi -pi]);
% % yticklabels({'-0.5\pi' 'pi'});
% % title('Phase')
% xlabel('\omega [rad/s]')
% 
% print -djpeg FRF.jpg -r300
%}

% % END plot of FRF
%% Least squares solution & Tikhonov solution
%
p_LSw = zeros(n_p,N/2);
p_TIKHw = zeros(n_p,N/2);
p_TIKHww = zeros(n_p,N/2);
L = eye(n_p); % zeroth order
% L = [1 -1]; % first order

for ind = 2:N2
    
    Aind = squeeze(A(:,:,ind));
    % Without regularization (least squares solution)
    p_LSw(:,ind) = ((Aind).'*Aind)\(Aind.')*dw(:,ind);  % left out eye(n_p)                      

% % LAMBDA FROM L CURVE
%
     [Utemp,stemp,~] = csvd(A(:,:,ind)); % for every frequency
     [reg_corner,rho,eta,reg_param] = l_curve(Utemp',stemp,S_p.'*pw(:,ind),'Tikh',L);
     lambda = reg_corner;
%}
%     lambda = (10000000/(freqRad(ind)^1.5))*trace((Aind.'*Aind)/trace(L.'*L));
%     lambda = 1;
%     lambda = lambda*trace((Aind.'*Aind)/trace(L.'*L)); 
      p_TIKHw(:,ind) = (Aind.'*Aind + lambda^2*(L.'*L))\Aind.'*dw(:,ind); 
      
      lambda = 100000/(freqRad(ind)^4)*trace((Aind.'*Aind)/trace(L.'*L));
      p_TIKHww(:,ind) = (Aind.'*Aind + lambda^2*(L.'*L))\Aind.'*dw(:,ind);
end

p_TIKHw = [p_TIKHw zeros(n_p,1) conj(p_TIKHw(:,end:-1:2))];
p_TIKH = real(ifft(p_TIKHw,[],2));
p_TIKHww = [p_TIKHww zeros(n_p,1) conj(p_TIKHww(:,end:-1:2))];
p_TIKHH = real(ifft(p_TIKHww,[],2));    % For use in question 4. 
p_LSw = [p_LSw zeros(n_p,1) conj(p_LSw(:,end:-1:2))];%
p_LS = real(ifft(p_LSw,[],2));
savefile = 'IdentifiedLoads.mat';
save(savefile, 'p_TIKHH');
% clear Utemp stemp
%}
%% SOLVE AND PLOT FOR VARYING LAMBDA
%{

L = eye(n_p); % zeroth order
kappa = zeros(N/2,1);
% L = [1 -1]; % first order
term = [1 2 4];


p_TIKHw = zeros(2,10000,length(term));
p_TIKH = zeros(2,10000,length(term));
p_LSw = zeros(n_p,N/2);
     
for j = 1:length(term)
    p_TIKHwt = zeros(n_p,N/2);
    
    for ind = 2:N2
        Aind = squeeze(A(:,:,ind));
        lambda = 100000/(freqRad(ind)^term(j))*trace((Aind.'*Aind)/trace(L.'*L)); %%%CHANGE!
        p_TIKHwt(:,ind) = (Aind.'*Aind + lambda^2*(L.'*L))\Aind.'*dw(:,ind);         
    end
    

p_TIKHw(:,:,j) = [p_TIKHwt zeros(n_p,1) conj(p_TIKHwt(:,end:-1:2))];
p_TIKH(:,:,j) = real(ifft(p_TIKHw(:,:,j),[],2));

end

figure('Name','LSsolution','units','normalized','outerposition',...
         [0 0 .5 .6]);
subplot(2,1,1)     
for j =1:length(term)
   plot(t,p_TIKH(1,:,j)); hold on;
end
% pp = S_p.'*p;
% plot(t,pp(2,:));
title('F_{wave}')
legend('10^5/(\omega^1)','10^5/(\omega^2)','10^5/(\omega^4)'); %'pw', 
xlabel('time [s]')
ylabel('force [N]')
hold off

subplot(2,1,2)     
for j =1:length(term)
   plot(t,p_TIKH(2,:,j)); hold on;
end
% pp = S_p.'*p;
% plot(t,pp(2,:));
title('F_{wind}')
% legend('10^5/(\omega^1)','10^5/(\omega^2)','10^5/(\omega^4)'); %'pw', 
xlabel('time [s]')
ylabel('force [N]')
hold off

 print -djpeg tikhwithlambdaparameter.jpg -r300

    
    
    
%}
% % Plot of the solutions
%{
% % End of force identification plot
figure('Name','Lcurvesolution','units','normalized','outerposition',...
         [0 0 .7 .6]);

for ind = 1:n_p
    
    subplot(2,1,ind); % Time domain
    Preal = (S_p(:,ind).'*p) - mean(S_p(:,ind).'*p);
    Plambda = (p_TIKHH(ind,:) - mean(p_TIKHH(ind,:)));
    
    plot(t,Preal,t,Plambda); %
    %ylim([-10^6;10^6]);
    xlim([0;10]);
    legend('p_{real}','p_{tikH\lambda*}');
    if ind == 1
        title('F_{wave}')
    else
        title('F_{wind}')
    end
    xlabel('time [s]')
    ylabel('force [N]')
        
%     subplot(2,2,ind+2); % Frequency domain
%     plot(freqRad(1:N2),S_p(:,ind).'*abs(pw(:,1:N2)),freqRad(1:N2),abs(p_LSw(ind,1:N2)),freqRad(1:N2),abs(p_TIKHw(ind,1:N2)),freqRad(1:N2),abs(p_TIKHww(ind,1:N2)))
%     legend('p_{real}','p_{LS}','p_{TikhLcurve}','p_{tikH\lambda}'); 
%   
%     xlim([1,freqRad(N/2)]);
%     xlabel('\omega [rad/s]')
%     ylabel('amplitude')

end

%  print -djpeg finalloadcomparisonzoom.jpg -r300
%}

% % plot of LS solution
%{

figure('Name','LSsolution','units','normalized','outerposition',...
         [0 0 .5 .5]);
     
subplot(2,1,1)
plot(t,p_TIKH(1,:),t,p_TIKH(2,:));
legend(['F_{wave,LS}';'F_{wind,LS}']);
ylim([-10^6;10^6]);
xlabel('time [s]')
ylabel('Force [N]')

subplot(2,1,2)
plot(freqRad(1:N2),abs(p_TIKHw(1,1:N2)),...
    freqRad(1:N2),abs(p_TIKHw(2,1:N2)));
ylabel('amplitude')
xlabel('\omega [rad/s]')
xlim([0 , freqRad(N2)]);
ylim([0 , 10^8]);
% 
% print -djpeg LSquares.jpg -r300

%}
