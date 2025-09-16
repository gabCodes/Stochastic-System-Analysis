%Time domain analysis, using example 7.1 as reference. 
ss_generator
dt = 0.05; fs = 1/dt;
T = 200; t = 0:dt:T; N = length(t);
nn = zeros(1,N);
w3 = randn(1,N)/sqrt(dt); %divide by sqrt(dt) due to lsim characteristics 

%Vertical turbulence only
u = [nn' nn' w3'];
%% Undamped system timeseries
y = lsim(A,B,C,D,u,t);
figure
sgtitle('Undamped system timeseries')
subplot(5, 1 , 1) ;
plot(t,y(:,1))
xlabel('time [s]'); ylabel('u / V [ -]') ; title('airspeed deviation');
subplot(5, 1 , 2 );
plot(t,y(:,2) * 180/pi)
xlabel(' time [s] ') ; ylabel ('\alpha [deg]') ; title ('angle of attack') ;
subplot(5, 1 , 3);
plot(t, y(:,3) * 180/pi)
xlabel('time [s]'); ylabel('\theta [deg]') ; title ('pitch angle');
subplot(5, 1 , 4);
plot(t, y(:,4) * 180/pi)
xlabel('time [s]') ; ylabel ('qc /V [deg]') ; title ('pitch rate');
subplot(5, 1 , 5);
plot(t, y(:,8))
xlabel('time [s]') ; ylabel ('n_z [-]') ; title ('load factor');
%% Damped system timeseries
y_pd = lsim(ss_ac_pd,u,t);
figure
sgtitle('Damped system timeseries')
subplot(5, 1 , 1) ;
plot(t,y_pd(:,1))
xlabel('time [s]'); ylabel('u / V [ -]') ; title('airspeed deviation');
subplot(5, 1 , 2 );
plot(t,y_pd(:,2) * 180/pi)
xlabel(' time [s] ') ; ylabel ('\alpha [deg]') ; title ('angle of attack') ;
subplot(5, 1 , 3);
plot(t, y_pd(:,3) * 180/pi)
xlabel('time [s]'); ylabel('\theta [deg]') ; title ('pitch angle');
subplot(5, 1 , 4);
plot(t, y_pd(:,4) * 180/pi)
xlabel('time [s]') ; ylabel ('qc /V [deg]') ; title ('pitch rate');
subplot(5, 1 , 5);
plot(t, y_pd(:,8))
xlabel('time [s]') ; ylabel ('n_z [-]') ; title ('load factor');


