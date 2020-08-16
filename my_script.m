
%

clear % clear all variables
%
%*********
%
ap= [   0.9759   3.5663    
         -0.0071  0.3810 ]

bp = [0.1699   0    
      0        0.0429]

cp = [    -4.2732   0
         0   0.0178 ]

dp = 0*ones(2,2)


%*********
%
% Natural Modes: Poles (Eigenvalues), Eigenvectors
%
[evec,eval] = eig(ap)   % evec contains eigenvectors
                        % eval contains poles or eigenvalues
                 
%******
%
%
t = [0:0.1:12];
u = [0*t' 0*t'];         % Set input u to zero for all time in order to generate zero input response;
                         % i.e. response to an initial condition x_o.

%
% Excite F404 (slow) temperature mode.
% This mode 
%            is associated with a pole at s = - 0.4.
%            has a time constant of  2.5 sec and settles in about 12.5 sec.
%            is the slowest of the F404's three (3) modes.
%            is associated with T_45
% Comment: It takes a long time to change temperature!
%
y = lsim(ss(ap, bp, eye(2,2), 0*ones(2,2)), u, t, evec(:,1)); 
plot(t,y)
grid
ylabel('States')
xlabel('Time (seconds)')
pause

%
y = lsim(ss(ap, bp, eye(2,2), 0*ones(2,2)), u, t, evec(:,2)); 
plot(t,y)
grid
title('F404 Slow RPM Mode: x_o = [ -0.9638 -0.2245 0.1440 ]')
ylabel('States')
xlabel('Time (seconds)')
pause

%

%
y = lsim(ss(ap, bp, eye(2,2), 0*ones(2,2)), u, t, evec(:,2)); 
plot(t,y)
grid
title('F404 Fast RPM Mode: x_o = [ 0.8676  -0.4781  -0.1368 ]')
ylabel('States')
xlabel('Time (seconds)')
pause


%return


%*********
%
% Transmission Zeros
%
%z = tzero(ss(ap,bp,cp,dp))                % transmission zeros
[p,z] = pzmap(ss(ap,bp,cp,dp))                % transmission zeros
%zdir =null([z*eye(2) - ap  -bp; cp dp])   % transmission zero directions

%*********
%
%
% F404 Engine TRANSFER FUNCTIONS: From u_i to x_j
%
sys = zpk(ss(ap,bp,eye(2,2),0*ones(2,2))) % Zeros, Poles, and Gains fron u_i to x_j



%*********
%
% Controllability 
%
cm = [bp ap*bp (ap^2)*bp]  % Controllability Matrix
rcm= rank(cm)              % Rank of Controllability Matrix


%*********
%
% Observability
%
om = [cp; 
      cp*ap;
      cp*(ap^2) ]          % Observability Matrix
rom = rank(om)             % Rank of Observability Matrix

%return

%*********
%
%
% F404 FREQUENCY RESPONSE: Unscaled Singular Values
%
% u = [ W_f (pph)  A_8 (sq in) ]
% x = [ N_2 (rpm)  N_25 (rpm)  T_45 (deg Fah) ]
% y = [ N_2 (rpm)              T_45 (deg Fah) ]
%
w = logspace(-2,3,100);
sv = sigma(ss(ap, bp, cp, dp),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Outputs: N_2, T_45; Inputs: W_f/110, A_8/22')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

%return

%*********
%
% F404 SVD ANALYSIS at DC
%
dc =  cp*inv(-ap)*bp
[udc,sdc,vdc] = svd(dc)


%*********
%
% F404 SVD Analysis at w = 0.1 rad/sec
%
%s = j*0.1
%g = cp*inv(s*eye(3)-ap)*bp + dp
%[u, s, v ] = svd(g)

%*********
%
% Step Response
%
t = [0:0.1:12];
y = step(ss(ap, bp, eye(2,2), 0*ones(2,2)), t); % Step response for each state and each input

%
% Step Response due to W_f
%

%
% N_2: Due to W_f
%
plot(t,y(:,1,1))
grid
title('N_2 response for W_f =  Unit Step')
ylabel('N_2 (rpm)')
xlabel('Time (seconds)')
pause
%return

%
% N_25: Due to W_f
%
plot(t,y(:,2,1))
grid
title('N_{25} response for W_f =  Unit Step')
ylabel('N_{25} (rpm)')
xlabel('Time (seconds)')
pause
%return

%
% T_45: Due to W_f
%
plot(t,y(:,2,1))
grid
title('T_{45} response for W_f =  Unit Step')
ylabel('T_{45} (deg F)')
xlabel('Time (seconds)')
pause
%return

%
% Step Response due to A_8
%

%
% N_2: Due to A_8
%
plot(t,y(:,1,2))
grid
title('N_2 response for A_8 =  Unit Step')
ylabel('N_2 (rpm)')
xlabel('Time (seconds)')
pause
%return

%
% N_25: Due to W_f
%
plot(t,y(:,2,2))
grid
title('N_{25} response for A_8 =  Unit Step')
ylabel('N_{25} (rpm)')
xlabel('Time (seconds)')
pause
%return

%
% T_45: Due to W_f
%
plot(t,y(:,2,2))
grid
title('T_{45} response for A_8 =  Unit Step')
ylabel('T_{45} (deg F)')
xlabel('Time (seconds)')
pause
%return

%return

%
% FACTS ON SCALING: Scaling affects the shape of singular value plots.
%                   It does not alter pole locations, zero locations.
%                   It does alter directionality information.
%

%
% Scaling Matrices
%
%  unew = su uold
%  xnew = sx xold
%  ynew = sy yold
%
su = diag( [1/110, 1/22] )
sx = diag( [1/250,  1/28] ) %1/350, 
sy = diag( [1/250,  1/28] )

%
% Scaled F404 Dynamics
%
%
% u = [ W_f/110    A_8/22 ]
% x = [ N_2/250  N_25/350  T_45/28 ]
% y = [ N_2/250            T_45/28 ]
%
ap = sx*ap*inv(sx)
bp = sx*bp*inv(su)
cp = sy*cp*inv(sx)
dp = sy*dp*inv(su)
    
%
% F404 FREQUENCY RESPONSE: Scaled Singular Values
%
sv = sigma(ss(ap, bp, cp, dp),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Outputs: N_2/250, T_45/350; Inputs: W_f/110, A_8/22')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause


%
% Augment Plant with Integrators at Plant Input and Plot Singular Values
%
[ns nc] = size(bp);                      % ns = number of inputs;  nc = number of controls;   
a = [ ap             bp
      0*ones(nc,ns)    0*ones(nc,nc) ]

b = [ 0*ones(ns,nc)
      eye(nc)      ]

c = [ cp  0*ones(nc,nc) ]

d = 0*ones(nc,nc)
sv = sigma(ss(a, b, c, d),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Design Plant Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

%return


%
%***********
%


%
% Design of Target Loop Singular Values Using Kalman Filter
%
ll =  inv(cp*inv(-ap)*bp + dp);     % Choose ll and lh to match singular values at all frequencies
lh = -inv(ap)*bp*ll;
l = [lh 
     ll];                           % ll, lh - for low and high frequency loop shaping

sv = sigma(ss(a, l, c, d),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Filter Open Loop (G_{FOL}) Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause



pnint = eye(nc)                                    % Process Noise Intensity Matrix
mu = 0.01;                                         % Measurement Noise Intesity; Used to adjust Kalman Filter Bandwidth
                                                   % Small mu - expensive sensor   - large bandwidth
                                                   % Large mu - inexpensive sensor - small bandwidth
mnint = mu*eye(nc)                                 % Measurement Noise Intensity Matrix 
% sysKal=ss(a, [b l], c, [d 0*ones(nc,nc)]);
% [kest, h, sig]= kalman(sysKal,pnint, mnint);  % Compute Filter Gain Matrix h
%[sig, poles, g1, rr] = care(a',c',l*l', mnint);                          
[sig, poles, g1] = care(a',c',l*l', mnint);                          
% Alternate Method for Computing h
h = g1';
sv = sigma(ss(a, h, c, d),w);
tsv = 20*log10(sv);
semilogx(w, tsv)
%clear sv
title('Target Loop (G_{KF}) Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

tolpoles = eig(a)                           % Target Open Loop Poles
%targzeros = tzero(a,h,c,0*ones(nc,nc))      % Target Open Loop Zeros
[targpoles,targzeros] = pzmap(ss(a,h,c,0*ones(nc,nc)))      % Target Open Loop Zeros
tclpoles = eig(a-h*c)                       % Target Closed Loop Poles

sv = sigma(ss(a-h*c, h, -c, eye(nc)),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Target Sensitivity (S_{KF}) Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause


sv = sigma(ss(a-h*c, h, c, 0*eye(nc)),w);
sv = 20*log10(sv);
semilogx(w, sv, w, 20*log10(10./w))
%clear sv
title('Target Complementary (T_{KF}) Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause


%
% Recover Target Loop By Solving Cheap LQR Problem
%
q = c'*c;                                            % State Weighting Matrix
rho = 1e-3;                                          % Cheap control recovery parameter;
                                                     % The smaller the parameter, the better the recovery.
r = rho*eye(nc)                                      % Control Weigthing Matrix
%[k, poles, g, rr] = care(a,b,q,r);                   % Compute Control Gain Matrix G
[k, poles, g] = care(a,b,q,r);                   % Compute Control Gain Matrix G


%
% Compensator Analysis
%
ak = [ a-b*g-h*c  0*ones(ns+nc,nc)
       g          0*ones(nc,nc) ]

bk = [ h
       0*ones(nc,nc) ]

ck = [0*ones(nc, ns+nc) eye(nc,nc) ]


%cpoles = eig(ak)                               % Compensator Poles
%czeros = tzero(a, h, g, 0*ones(nc,nc))         % Compensator Zeros
[cpoles, czeros] = pzmap(ss(a, h, g, 0*ones(nc,nc)))         % Compensator Zeros
%zerocheck = tzero(ak, bk, ck, 0*ones(nc,nc))   % Check Compensator Zeros
[polecheck, zerocheck] = pzmap(ss(ak, bk, ck, 0*ones(nc,nc)))   % Check Compensator Zeros

sv = sigma(ss(ak, bk, ck, 0*eye(nc)),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Compensator Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause


%
% Open Loop Analysis
%
al = [ ap                     bp*ck
       0*ones(ns+nc+nc,ns)    ak    ]

bl = [ 0*ones(ns,nc)
       bk ]
    
cl = [ cp  0*ones(nc,ns+nc+nc) ]
    
%olpoles = eig(al)                          % Open Loop Poles
%olzeros = tzero(al,bl,cl,0*ones(nc,nc))    % Open Loop Zeros
[olpoles, olzeros] = pzmap(ss(al,bl,cl,0*ones(nc,nc)))    % Open Loop Zeros    
sv = sigma(ss(al, bl, cl, 0*eye(nc)),w);
sv = 20*log10(sv);
semilogx(w, sv, w, tsv)
%clear sv
title('Open Loop Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause   

%
% Closed Loop Analysis
%
clpoles = eig(al-bl*cl)           % Closed Loop Poles
clpkf = eig(a - h*c)              % Closed Loop Poles Due to Kalman Filter
clpreg = eig(a - b*g)             % Closed Loop Poles Due to Regulator


sv = sigma(ss(al-bl*cl, bl, -cl, eye(nc)),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Sensitivity Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause 

sv = sigma(ss(al-bl*cl, bl, cl, 0*eye(nc)),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Complementary Sensitivity Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause 

%
% Step Response in Closed Loop 
%

[y,t] = step(ss(al-bl*cl, bl, cl, 0*eye(nc)));
plot(t,y(:,1,1))
xlabel('Time (s)')
ylabel('Amplitude')
title('output 1 response caused by input 1')
pause
%return

%
plot(t,y(:,1,2))
xlabel('Time (s)')
ylabel('Amplitude')
title('output 1 response caused by input 2')
pause
%return

plot(t,y(:,2,1))
xlabel('Time (s)')
ylabel('Amplitude')
title('output 2 response caused by input 1')
pause
%return

%
plot(t,y(:,2,2))
xlabel('Time (s)')
ylabel('Amplitude')
title('output 2 response caused by input 2')
pause
%return
disp('You can see a good tracking and a good disturbance rejection')

