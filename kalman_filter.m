function [x_filter,P_filter,diff,Kmagn] = kalman_filter(A,B,C,D,Q,R,P0,x0,p,d);

n_s = size(A,1);	% Number of states
n_d = size(d,1);	% Number of observations
N = size(d,2);		% Number of samples
x_filter = zeros(n_s,N);		
P_filter = zeros(n_s,n_s,N);
diff = zeros(n_d,N);
Kmagn = zeros(N,1);
for time = 1:N
	%	Time update
	if (time==1)
		x_filter(:,time) = A*x0; 												
		P_filter(:,:,time) = A*P0*A' + Q;
	else
		x_filter(:,time) = A*x_filter(:,time-1)+B*p(:,time-1); 												
		P_filter(:,:,time) = A*P_filter(:,:,time-1)*A' + Q;
	end
	
	%	Measurement update
	K_k = P_filter(:,:,time)*C'*inv(C*P_filter(:,:,time)*C'+R);
    x_filter(:,time) = x_filter(:,time) + K_k*(d(:,time)-C*x_filter(:,time)-D*p(:,time));
    diff(:,time) = d(:,time)-C*x_filter(:,time)-D*p(:,time);
	P_filter(:,:,time) = P_filter(:,:,time) - K_k*C*P_filter(:,:,time);
    Kmagn(time) = norm(K_k);
end