clear all;
clc;
rng(1)
% State system generator based off of cit2s, outputs the state space model
% corresponding to the turbulence intensity & scale length.

% INPUT TURBULENCE- AND AIRCRAFT PARAMETERS

% AIRCRAFT FLIGHT CONDITION 'LANDING'.
V     = 51.4;
m     = 4556;
twmuc = 2*76;
KY2   = 0.980;
c     = 2.022;
S     = 24.2;
lh    = 5.5;
g     = 9.80665;

% TURBULENCE PARAMETERS
sigma = 2;
Lg    = 1500;

sigmaug_V = sigma/V;
sigmaag   = sigma/V;

% AIRCRAFT SYMMETRIC AERODYNAMIC DERIVATIVES : 
CX0 = 0.0000;     CZ0  =-1.1360;     Cm0  =  0.0000;
CXu =-0.2173;     CZu  =-2.2720;     Cmu  =  0.0000;
CXa = 0.4692;     CZa  =-5.1600;     Cma  = -0.4000;
CXq = 0.0000;     CZq  =-3.8400;     Cmq  = -7.3500;
CXd = 0.0000;     CZd  =-0.6238;     Cmd  = -1.5530;
CXfa= 0.0000;     CZfa =-1.4050;     Cmfa = -3.6150;
                  CZfug= 0.0000;     Cmfug= -Cm0*lh/c;
                  CZfag= CZfa-CZq;   Cmfag=  Cmfa-Cmq;

% CALCULATION OF AIRCRAFT SYMMETRIC STABILITY DERIVATIVES
xu   = (V/c)*(CXu/twmuc);
xa   = (V/c)*(CXa/twmuc);
xt   = (V/c)*(CZ0/twmuc);
xq   = 0;
xd   = (V/c)*(CXd/twmuc);
xug  = xu;
xfug = 0;
xag  = xa;
xfag = 0;

zu   = (V/c)*( CZu/(twmuc-CZfa));
za   = (V/c)*( CZa/(twmuc-CZfa));
zt   = (V/c)*(-CX0/(twmuc-CZfa));
zq   = (V/c)*((CZq+twmuc)/(twmuc-CZfa));
zd   = (V/c)*( CZd/(twmuc-CZfa));
zug  = zu;
zfug = (V/c)*( CZfug/(twmuc-CZfa));
zag  = za;
zfag = (V/c)*( CZfag/(twmuc-CZfa));

mu   = (V/c)*(( Cmu+CZu*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
ma   = (V/c)*(( Cma+CZa*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
mt   = (V/c)*((-CX0*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
mq   = (V/c)*(Cmq+Cmfa*(twmuc+CZq)/(twmuc-CZfa))/(twmuc*KY2);
md   = (V/c)*((Cmd+CZd*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
mug  = mu;
mfug = (V/c)*(Cmfug+CZfug*Cmfa/(twmuc-CZfa))/(twmuc*KY2);
mag  = ma;
mfag = (V/c)*(Cmfag+CZfag*Cmfa/(twmuc-CZfa))/(twmuc*KY2);

% STATE- AND INPUT MATRICES
A=[xu xa xt 0    xug                  xag       0;
   zu za zt zq   zug-zfug*V/Lg*(c/V)  zag       zfag*(c/V);
   0  0  0  V/c  0                    0         0;
   mu ma mt mq   mug-mfug*V/Lg*(c/V)  mag       mfag*(c/V);
   0  0  0  0   -V/Lg                 0         0;
   0  0  0  0    0                    0         1;
   0  0  0  0    0                   -(V/Lg)^2 -2*V/Lg];

B=...
 [xd 0                                 0;
  zd zfug*(c/V)*sigmaug_V*sqrt(2*V/Lg) zfag*(c/V)*sigmaag*sqrt(3*V/Lg);
  0  0                                 0;
  md mfug*(c/V)*sigmaug_V*sqrt(2*V/Lg) mfag*(c/V)*sigmaag*sqrt(3*V/Lg);
  0  sigmaug_V*sqrt(2*V/Lg)            0;
  0  0                                 sigmaag*sqrt(3*V/Lg);
  0  0                                 (1-2*sqrt(3))*sigmaag*sqrt((V/Lg)^3)];

C=...
    [1     0     0     0     0     0     0;
     0     1     0     0     0     0     0;
     0     0     1     0     0     0     0;
     0     0     0     1     0     0     0;
     0     0     0     0     1     0     0;
     0     0     0     0     0     1     0;
     0     0     0     0     0     0     1];
D=...
    [0     0     0;
     0     0     0;
     0     0     0;
     0     0     0;
     0     0     0;
     0     0     0;
     0     0     0];

ss_ac = ss(A,B,C,D);

%% Extending the output equation to include the load factor
%Using definition of n(t) from Example 7.3
C_nz = (-V/g)*(ss_ac.A(2,:)-ss_ac.A(3,:));
D_nz = (-V/g)*(ss_ac.B(2,:)-ss_ac.B(3,:));
C = [C; C_nz];
D = [D; D_nz];

ss_ac = ss(A,B,C,D);
%Ktheta val = 0.11637, negative feedback
ss_pd = ss(A,-B(:,1),[0 0 1 0 0 0 0], 0);
%sisotool(ss_pd)
Kt = -0.11637;
K = [0 0 Kt 0 0 0 0];
A_pd = A-B(:,1)*K;
ss_ac_pd = ss(A_pd, B, C, D);
disp('State Spaces has been loaded');
