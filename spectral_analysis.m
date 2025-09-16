time_domain_analysis
%% Analytical PSD Calculation
% COMPUTE FREQUENCY RESPONSE FUNCTION AND PSD
w = logspace(-2,2,4000);
mag = bode(A,B,C(1,:),D(1,:),3,w); Suu_u   = mag.*mag;
mag = bode(A,B,C(2,:),D(2,:),3,w); Saa_u   = mag.*mag;
mag = bode(A,B,C(3,:),D(3,:),3,w); Stt_u   = mag.*mag;
mag = bode(A,B,C(4,:),D(4,:),3,w); Sqq_u   = mag.*mag;
mag = bode(A,B,C(8,:),D(8,:),3,w); Snz_u   = mag.*mag;

mag = bode(A_pd,B,C(1,:),D(1,:),3,w); Suu_d   = mag.*mag;
mag = bode(A_pd,B,C(2,:),D(2,:),3,w); Saa_d   = mag.*mag;
mag = bode(A_pd,B,C(3,:),D(3,:),3,w); Stt_d   = mag.*mag;
mag = bode(A_pd,B,C(4,:),D(4,:),3,w); Sqq_d   = mag.*mag;
mag = bode(A_pd,B,C(8,:),D(8,:),3,w); Snz_d   = mag.*mag;

Sa_u = [Suu_u Saa_u Stt_u Sqq_u Snz_u];
Sa_d = [Suu_d Saa_d Stt_d Sqq_d Snz_d];

clf
subplot(3,2,1) 
loglog(w,Sa_u(:,1),'-r', w, Sa_d(:,1),'-b', 'LineWidth',1);
axis(10.^[-2,2,-15,0])
xlabel('\omega [rad/s]'); ylabel('Suu [1/(rad/s)]')
grid on
legend('Undamped', 'Damped')

subplot(3,2,2)
loglog(w,Sa_u(:,2),'-r', w, Sa_d(:,2), '-b','LineWidth',1);
axis(10.^[-2,2,-15,0])
xlabel('\omega [rad/s]'); ylabel('Saa [rad^2/(rad/s)]')
grid on
legend('Undamped', 'Damped')

subplot(3,2,3)
loglog(w,Sa_u(:,3),'-r', w, Sa_d(:,3), '-b','LineWidth',1);
axis(10.^[-2,2,-15,0])
xlabel('\omega [rad/s]'); ylabel('Stt [rad^2/(rad/s)]')
grid on
legend('Undamped', 'Damped')

subplot(3,2,4)
loglog(w,Sa_u(:,4),'-r', w, Sa_d(:,4), '-b','LineWidth',1);
axis(10.^[-2,2,-20,-5])
xlabel('\omega [rad/s]'); ylabel('Sqq [1/(rad/s)]')
grid on
legend('Undamped', 'Damped')

subplot(3,2,5)
loglog(w,Sa_u(:,5),'-r', w, Sa_d(:,5), '-b','LineWidth',1);
axis(10.^[-2,2,-10,0])
xlabel('\omega [rad/s]'); ylabel('Snznz [1/rad/s]')
grid on
legend('Undamped', 'Damped')
disp('Analytical PSD plotted, press any key to proceed with the FFT ones')
pause

%% Experimental PSD (FFT)
% DEFINE FREQUENCY VECTOR FOR PLOTTING
omega = 2*pi*(fs/N)*(0:(N/2)-1);

%Selecting the relevant states
hatu   = y(:,1);    hatu_pd     = y_pd(:,1);
alpha  = y(:,2);    alpha_pd    = y_pd(:,2);
theta  = y(:,3);    theta_pd    = y_pd(:,3);
qcV    = y(:,4);    qcV_pd      = y_pd(:,4);
Nz     = y(:,8);    Nz_pd       = y_pd(:,8);


% FFT ALL SIGNALS
U      = dt*fft(hatu);      U_pd     =  dt*fft(hatu_pd);
ALPHA  = dt*fft(alpha);     ALPHA_pd =  dt*fft(alpha_pd);  
THETA  = dt*fft(theta);     THETA_pd =  dt*fft(theta_pd); 
QCV    = dt*fft(qcV);       QCV_pd   =  dt*fft(qcV_pd);
NZ     = dt*fft(Nz);        NZ_pd    =  dt*fft(Nz_pd);

% COMPUTE PSDs, note we take half the vectors to obtain only the positive
% frequencies (the remaining values correspond to Nyquist + the negative frequencies
Pu      = real((1/T)*     U.*conj(U));      Pu      = Pu(1:round(N/2)-1);
Palpha  = real((1/T)* ALPHA.*conj(ALPHA));  Palpha  = Palpha(1:round(N/2)-1);
Ptheta  = real((1/T)* THETA.*conj(THETA));  Ptheta = Ptheta(1:round(N/2)-1);
PqcV    = real((1/T)*   QCV.*conj(QCV));    PqcV = PqcV(1:round(N/2)-1);
Pnz     = real((1/T)*    NZ.*conj(NZ));     Pnz = Pnz(1:round(N/2)-1);

Pu_pd      = real((1/T)*     U_pd.*conj(U_pd));     Pu_pd     = Pu_pd(1:round(N/2)-1);
Palpha_pd  = real((1/T)* ALPHA_pd.*conj(ALPHA_pd)); Palpha_pd = Palpha_pd(1:round(N/2)-1);
Ptheta_pd  = real((1/T)* THETA_pd.*conj(THETA_pd)); Ptheta_pd = Ptheta_pd(1:round(N/2)-1);
PqcV_pd    = real((1/T)*   QCV_pd.*conj(QCV_pd));   PqcV_pd   = PqcV_pd(1:round(N/2)-1);
Pnz_pd     = real((1/T)*    NZ_pd.*conj(NZ_pd));    Pnz_pd    = Pnz_pd(1:round(N/2)-1);

Sf_u = [Pu, Palpha, Ptheta, PqcV, Pnz];
Sf_d = [Pu_pd, Palpha_pd, Ptheta_pd, PqcV_pd, Pnz_pd];

clf
subplot(3,2,1) 
loglog(omega,Sf_u(:,1),'-r', omega, Sf_d(:,1), '-b', 'LineWidth',1);
axis(10.^[-2,2,-15,0])
xlabel('\omega [rad/s]'); ylabel('Suu [1/(rad/s)]')
grid on
legend('Undamped', 'Damped')

subplot(3,2,2)
loglog(omega,Sf_u(:,2),'-r', omega, Sf_d(:,2), '-b','LineWidth',1);
axis(10.^[-2,2,-15,0])
xlabel('\omega [rad/s]'); ylabel('Saa [rad^2/(rad/s)]')
grid on
legend('Undamped', 'Damped')

subplot(3,2,3)
loglog(omega,Sf_u(:,3),'-r', omega, Sf_d(:,3), '-b','LineWidth',1);
axis(10.^[-2,2,-15,0])
xlabel('\omega [rad/s]'); ylabel('Stt [rad^2/(rad/s)]')
grid on
legend('Undamped', 'Damped')

subplot(3,2,4)
loglog(omega,Sf_u(:,4),'-r', omega, Sf_d(:,4), '-b','LineWidth',1);
axis(10.^[-2,2,-20,-5])
xlabel('\omega [rad/s]'); ylabel('Sqq [1/(rad/s)]')
grid on
legend('Undamped', 'Damped')

subplot(3,2,5)
loglog(omega,Sf_u(:,5),'-r', omega, Sf_d(:,5), '-b','LineWidth',1);
axis(10.^[-2,2,-10,0])
xlabel('\omega [rad/s]'); ylabel('Snznz [1/rad/s]')
grid on
legend('Undamped', 'Damped')
disp('FFT PSD plotted, press any key to proceed with the Pwelch ones')
pause
%% Experimental PSD (pwelch)
% DEFINE FREQUENCY VECTOR FOR PLOTTING
omega = 2*pi*fs*(1:(N/2))/N;

% Selecting the relevant states
hatu   = y(:,1);    hatu_pd     = y_pd(:,1);
alpha  = y(:,2);    alpha_pd    = y_pd(:,2);
theta  = y(:,3);    theta_pd    = y_pd(:,3);
qcV    = y(:,4);    qcV_pd      = y_pd(:,4);
Nz     = y(:,8);    Nz_pd       = y_pd(:,8);

% Undamped
[Suu_pw fw] = pwelch(hatu,100,50,N,fs,'onesided');
[Saa_pw fw] = pwelch(alpha,100,50,N,fs,'onesided');
[Stt_pw fw] = pwelch(theta,100,50,N,fs,'onesided');
[Sqq_pw fw] = pwelch(qcV,100,50,N,fs,'onesided');
[Snz_pw fw] = pwelch(Nz,100,50,N,fs,'onesided');

Suu_pw = Suu_pw/2;          % Adjust for negative frequencies
Suu_pw = Suu_pw(2:round(N/2));   % Remove the zero frequency to align
Saa_pw = Saa_pw/2;          % Adjust for negative frequencies
Saa_pw = Saa_pw(2:round(N/2));   % Remove the zero frequency to align
Stt_pw = Stt_pw/2;          % Adjust for negative frequencies
Stt_pw = Stt_pw(2:round(N/2));   % Remove the zero frequency to align
Sqq_pw = Sqq_pw/2;          % Adjust for negative frequencies
Sqq_pw = Sqq_pw(2:round(N/2));   % Remove the zero frequency to align
Snz_pw = Snz_pw/2;          % Adjust for negative frequencies
Snz_pw = Snz_pw(2:round(N/2));   % Remove the zero frequency to align

% Damped
[Suu_pw_pd fw] = pwelch(hatu_pd,100,50,N,fs,'onesided');
[Saa_pw_pd fw] = pwelch(alpha_pd,100,50,N,fs,'onesided');
[Stt_pw_pd fw] = pwelch(theta_pd,100,50,N,fs,'onesided');
[Sqq_pw_pd fw] = pwelch(qcV_pd,100,50,N,fs,'onesided');
[Snz_pw_pd fw] = pwelch(Nz_pd,100,50,N,fs,'onesided');

Suu_pw_pd = Suu_pw_pd/2;          % Adjust for negative frequencies
Suu_pw_pd = Suu_pw_pd(2:round(N/2));   % Remove the zero frequency to align
Saa_pw_pd = Saa_pw_pd/2;          % Adjust for negative frequencies
Saa_pw_pd = Saa_pw_pd(2:round(N/2));   % Remove the zero frequency to align
Stt_pw_pd = Stt_pw_pd/2;          % Adjust for negative frequencies
Stt_pw_pd = Stt_pw_pd(2:round(N/2));   % Remove the zero frequency to align
Sqq_pw_pd = Sqq_pw_pd/2;          % Adjust for negative frequencies
Sqq_pw_pd = Sqq_pw_pd(2:round(N/2));   % Remove the zero frequency to align
Snz_pw_pd = Snz_pw_pd/2;          % Adjust for negative frequencies
Snz_pw_pd = Snz_pw_pd(2:round(N/2));   % Remove the zero frequency to align

Sp_u = [Suu_pw, Saa_pw, Stt_pw, Sqq_pw, Snz_pw];
Sp_d = [Suu_pw_pd, Saa_pw_pd, Stt_pw_pd, Sqq_pw_pd, Snz_pw_pd];

clf
subplot(3,2,1) 
loglog(omega,Sp_u(:,1),'-r', omega, Sp_d(:,1), '-b', 'LineWidth',1);
axis(10.^[-2,2,-15,0])
xlabel('\omega [rad/s]'); ylabel('Suu [1/(rad/s)]')
grid on
legend('Undamped', 'Damped')

subplot(3,2,2)
loglog(omega,Sp_u(:,2),'-r', omega, Sp_d(:,2), '-b','LineWidth',1);
axis(10.^[-2,2,-15,0])
xlabel('\omega [rad/s]'); ylabel('Saa [rad^2/(rad/s)]')
grid on
legend('Undamped', 'Damped')

subplot(3,2,3)
loglog(omega,Sp_u(:,3),'-r', omega, Sp_d(:,3), '-b','LineWidth',1);
axis(10.^[-2,2,-15,0])
xlabel('\omega [rad/s]'); ylabel('Stt [rad^2/(rad/s)]')
grid on
legend('Undamped', 'Damped')

subplot(3,2,4)
loglog(omega,Sp_u(:,4),'-r', omega, Sp_d(:,4), '-b','LineWidth',1);
axis(10.^[-2,2,-20,-5])
xlabel('\omega [rad/s]'); ylabel('Sqq [1/(rad/s)]')
grid on
legend('Undamped', 'Damped')

subplot(3,2,5)
loglog(omega,Sp_u(:,5),'-r', omega, Sp_d(:,5), '-b','LineWidth',1);
axis(10.^[-2,2,-10,0])
xlabel('\omega [rad/s]'); ylabel('Snznz [1/rad/s]')
grid on
legend('Undamped', 'Damped')

disp('PSDs Complete!')

%% Variances

%% Computed using timetraces and Var.m
var_u = var(y);
colNames = {'u/V [-]', '\alpha [rad^2]', 'theta [rad^2]', 'qc/V [rad^2]', 'NA', 'NA1', 'NA2', 'n_z [-]'};
% Type tab_u in console for Var undamped variances
tab_u = array2table(var_u, 'VariableNames',colNames);
var_pd = var(y_pd);
% Type tab_pd in console for Var damped variances
tab_pd = array2table(var_pd, 'VariableNames',colNames);

%% Computed using analytic PSD & basic crude integration scheme
colNames = {'u/V [-]', '\alpha [rad^2]', 'theta [rad^2]', 'qc/V [rad^2]', 'n_z [-]'};
var_au = zeros(1,5);
for j =1:5
    for i=1:length(w) - 1
        var_au(j) = var_au(j)+(w(i+1) - w(i))*Sa_u(i,j);
    end
end
var_au = var_au/pi;
% Type tab_au in console for analytic PSD undamped variances
tab_au = array2table(var_au, 'VariableNames',colNames);

var_ad = zeros(1,5);
for j =1:5
    for i=1:length(w) - 1
        var_ad(j) = var_ad(j)+(w(i+1) - w(i))*Sa_d(i,j);
    end
end
var_ad = var_ad/pi;
% Type tab_ad in console for analytic PSD damped variances
tab_ad = array2table(var_ad, 'VariableNames',colNames);

%% Computed using PSD estimates & basic crude integration scheme
% FFT estimates
var_fu = zeros(1,5);
for j =1:5
    for i=1:length(omega) - 1
        var_fu(j) = var_fu(j)+(omega(i+1) - omega(i))*Sf_u(i,j);
    end
end
var_fu = var_fu/pi;
% Type tab_fu in console for FFT PSD undamped variances
tab_fu = array2table(var_fu, 'VariableNames',colNames);

var_fd = zeros(1,5);
for j =1:5
    for i=1:length(omega) - 1
        var_fd(j) = var_fd(j)+(omega(i+1) - omega(i))*Sf_d(i,j);
    end
end
var_fd = var_fd/pi;
% Type tab_fd in console for FFT PSD damped variances
tab_fd = array2table(var_fd, 'VariableNames',colNames);

% Pwelch estimates
var_pwu = zeros(1,5);
for j =1:5
    for i=1:length(omega) - 1
        var_pwu(j) = var_pwu(j)+(omega(i+1) - omega(i))*Sp_u(i,j);
    end
end
var_pwu = var_pwu/pi;
% Type tab_pwu in console for pwelch PSD undamped variances
tab_pwu = array2table(var_pwu, 'VariableNames',colNames);

var_pwd = zeros(1,5);
for j =1:5
    for i=1:length(omega) - 1
        var_pwd(j) = var_pwd(j)+(omega(i+1) - omega(i))*Sp_d(i,j);
    end
end
var_pwd = var_pwd/pi;
% Type tab_pwd in console for FFT PSD damped variances
tab_pwd = array2table(var_pwd, 'VariableNames',colNames);