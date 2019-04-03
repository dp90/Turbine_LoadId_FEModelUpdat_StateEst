% % 
% % Model updating through gradient descent. Loss is defined as difference
% % between eigenfrequencies computed from FE model and eigenfrequencies
% % derived from acceleration data (DetEigFreqFromAccData.m). Change 2 input
% % parameters at a time for clarity. 
% % Goal update FE model for improved predictions on mechanical and 
% % dynamical behavior. 
% % 
% % -----------------------------------------------------------------------
% % -----------------------------------------------------------------------

% eigenvectors;

n_steps = 21;

pct_below = 0.5;
pct_above = 0.5;

Esoil_orig = 3e9;           % 1.744e9 (icm 6), 1.75e9 (icm 2)
Est_orig = 210e9;           % 202.8e9 (icm 6), 208.8 (icm 1)
rhost_orig = 7850;          % Optimum bij 8137 (icm 3,4,5), 7876 (icm 1), 7845 (icm 2)
Mtop_orig = 350000;         %

for i = 1:4
    for j = 1:4
        if i ~= j
            var1 = i;   % 1,2,3,4 for Esoil, Est, rhost, Mtop
            var2 = j;   % 1,2,3,4 for Esoil, Est, rhost, Mtop

            % Combination 2&3 is interesting. 

            switch var1
                case 1
                    x1 = linspace(Esoil_orig*(1-pct_below),Esoil_orig*(1+pct_above),n_steps);
                case 2
                    x1 = linspace(Est_orig*(1-pct_below),Est_orig*(1+pct_above),n_steps);
                case 3
                    x1 = linspace(rhost_orig*(1-pct_below),rhost_orig*(1+pct_above),n_steps);
                case 4
                    x1 = linspace(Mtop_orig*(1-pct_below),Mtop_orig*(1+pct_above),n_steps);
            end
            switch var2
                case 1
                    x2 = linspace(Esoil_orig*(1-pct_below),Esoil_orig*(1+pct_above),n_steps);
                case 2
                    x2 = linspace(Est_orig*(1-pct_below),Est_orig*(1+pct_above),n_steps);
                case 3
                    x2 = linspace(rhost_orig*(1-pct_below),rhost_orig*(1+pct_above),n_steps);
                case 4
                    x2 = linspace(Mtop_orig*(1-pct_below),Mtop_orig*(1+pct_above),n_steps);
            end

            assignin('base','x1',x1);
            assignin('base','x2',x2);

            RunOptim;
            costnew(i,j) = history.fval(end);

            sclx1 = evalin('base', 'sclx1');
            sclx2 = evalin('base', 'sclx2');
            x1old = x1(11)./sclx1;
            x2old = x2(11)./sclx2;
            f_old(i,j) = ObjFun( [x1old,x2old] );
            
            f_new = ObjFun( [history.x(end,1) , history.x(end,2)] );
            % improvement = f_old/f_new - 1
            % history.x(end,:)
            Omega_new(i,j,:) = evalin('base','Omega1');
            
        end
    end
end



Omega_meas = [0.34 0.86 5.04 6.64 9.83 12.62];
Omega_12 = [ Omega_new(1,2,1) Omega_new(1,2,2) Omega_new(1,2,3) Omega_new(1,2,5) Omega_new(1,2,6) Omega_new(1,2,10) ] / (2*pi);
Omega_13 = [ Omega_new(1,3,1) Omega_new(1,3,2) Omega_new(1,3,3) Omega_new(1,3,5) Omega_new(1,3,6) Omega_new(1,3,10) ] / (2*pi);
Omega_14 = [ Omega_new(1,4,1) Omega_new(1,4,2) Omega_new(1,4,3) Omega_new(1,4,5) Omega_new(1,4,6) Omega_new(1,4,10) ] / (2*pi);
Omega_23 = [ Omega_new(2,3,1) Omega_new(2,3,2) Omega_new(2,3,3) Omega_new(2,3,5) Omega_new(2,3,6) Omega_new(2,3,10) ] / (2*pi);
Omega_24 = [ Omega_new(2,4,1) Omega_new(2,4,2) Omega_new(2,4,3) Omega_new(2,4,5) Omega_new(2,4,6) Omega_new(2,4,10) ] / (2*pi);
Omega_34 = [ Omega_new(3,4,1) Omega_new(3,4,2) Omega_new(3,4,3) Omega_new(3,4,5) Omega_new(3,4,6) Omega_new(3,4,10) ] / (2*pi);

% 1 & 2: 0.0661
% 1 & 6: 0.0659
% 2 & 6: 0.0082
