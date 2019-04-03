% %             
% % State estimation with Kalman filter
% %
% % -----------------------------------------------------------------------
% % -----------------------------------------------------------------------

close all; clear; clc; 
FE_modelUPDATE;

load('Forces.mat');
load('IdentifiedLoads.mat');
load('Accelerations_B.mat');

reduced_measurements = 0;    % 0 is all measurements are included, 1 means reduced number of measurements
no_meas_incl = 2;            % number of measurements to include
if reduced_measurements == 1
    S_d = S_d(1:no_meas_incl,:);
    d = d(1:no_meas_incl,:);
end
dt = t(2) - t(1);
n_d = size(S_d,1);

p = [p_TIKHH];
c15=selectdof(DOF,15.01);   % Fwave
c35=selectdof(DOF,35.01);   % Fwind
S_p = [c15' c35'];          % [np x nDOF]

xi = 0.02;                       % Constant damping ratio assumed. Could be proportional
nDOF = length(Omega);

% State-space matrices A, B, C and D
S_v = zeros(size(S_d));          % selection matrix for velocity
S_a = S_d;                       % selection matrix for acceleration
S_d = zeros(size(S_v));
% ASSEMBLE SYSTEM MATRIX A
A_cm = [zeros(nDOF,nDOF) , eye(nDOF) ; -diag(Omega.^2) , -diag(2*xi*Omega)];
A = expm(A_cm*dt);
% ASSEMBLE INPUT MATRIX B
B_cm = [zeros(nDOF,size(S_p,2)); Phi'*S_p];
B = (A - eye(2*nDOF)) * inv(A_cm) * B_cm;	
% ASSEMBLE OUTPUT MATRIX C
C = [S_d*Phi - S_a*Phi*diag(Omega.^2) , S_v*Phi - S_a*Phi*diag(2*xi*Omega)];
% ASSEMBLE DIRECT TRANSMISSION MATRIX D
D = S_a*Phi*Phi'*S_p;

% State estimation
% w = 1e-9*ones(2*nDOF,1)                 % Process noise covariance matrix
% v = 1 * var(d,[],2);                
% R = v*v';                              % Measurement noise covariance matrix (estimate)

Transform = [Phi zeros(nDOF,nDOF) ; zeros(nDOF,nDOF) Phi];
Q = 0.01*eye(2*nDOF);                   % Process noise covariance matrix
R = 0.01*eye(n_d);                      % Measurement noise covariance matrix (estimate)             
x0 = 1*ones(2*nDOF,1);                  % Initial state estimate
P0 = 0.01*eye(2*nDOF);                  % Initial state error covariance matrix

[xm_filter,P_filter,diff,Kmagn] = kalman_filter(A,B,C,D,Q,R,P0,x0,p,d);
x_filter = Transform*xm_filter;

figure;
plot(t,xm_filter(1:2,:))
xlabel('Time [s]','fontsize',16);
ylabel('Modal coordinate z [-]','fontsize',16);
legend('First mode', 'Second mode');
print -djpeg ModalCoords.jpg -r300


%%
for i = 1:10000
    xnorm = norm(xm_filter(:,i));
    for j = 1:2
        xm_normalized(j,i) = xm_filter(j,i) / xnorm;
    end
end

NormDiff = zeros(1,10000);
for i = 1:10000
    NormDiff(i) = norm(diff(:,i));
end

% figure;
% plot(t,NormDiff)
% xlabel('Time [s]','fontsize',16);
% ylabel('Norm of difference between measurement and predicted measurement','fontsize',16);
% print -djpeg DiffNorm.jpg -r300
% 
% figure;
% plot(t,xm_normalized);
% xlabel('Time [s]','fontsize',16);
% ylabel('Modal coordinate z [-]','fontsize',16);
% legend('First mode', 'Second mode');
% print -djpeg ModalCoordsNorm.jpg -r300
% 
% figure;
% plot(t,x_filter(103,:));
% xlabel('Time [s]');
% ylabel('Displacement x [m]');
% legend('degree of freedom 103');

%%
% figure;
% animdisp(Nodes,Elements,Types,DOF,x_filter(1:105,:));

% plot(t,ci*x_filter,t,cii*x_filter,t,ciii*x_filter);
% title('displacement ');
% % axis([0 0.1 -50 200]);
% legend('displacment at tip', 'displacement at node 5','displacement at node 10','displacement at node 15');
% xlabel('Time [s]');
% ylabel('displacment [mm]');
% 
% clear ci cii ciii c19 c20 c21








