ss_generator
%% Undapmed aircraft characteristics
eig_ac = eig(ss_ac);
if all(real(eig_ac') < 0)
    disp('Undamped aircraft is stable');
else
    disp('Undamped aircraft is NOT stable');
end
%w_n_p = abs(eig_ac(3))  %natural frequecny
%zeta_p = - real(eig_ac(3))/w_n_p % damping coefficient zeta
%pzmap(ss_ac)

%% Damped aircraft characteristics
eig_ac_pd = eig(ss_ac_pd);
if all(real(eig_ac_pd') < 0)
    disp('Damped aircraft is stable');
else 
    disp('Damped aircraft is NOT stable');
end
%w_n_pd = abs(eig_ac_pd(3))  %natural frequecny
%zeta_pd = - real(eig_ac_pd(3))/w_n_pd % damping coefficient zeta
%pzmap(ss_ac_pd)
