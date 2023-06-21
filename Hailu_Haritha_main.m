 global stress_h;
 global Kg_sh;
 global shear_h;
 global zeta_c;

rho = 1050.0;  % kg/m^3
PI = 3.14159265;
Tol = 0.00000000001;

ro_h = 0.0003282;  % outer radius
h_h = 0.000021175;
Gc_h = 1.07;
Gm_h = 1.25;
Ge_1 = 2.19;
Ge_2 = 1.64;
stress_h = 120000.0;  % Pa
shear_h = 1.5;  % Pa
phi_e0 = 0.06;
phi_m0 = 0.09;
phi_f = 0.7;
phi_k0 = [0.11137, 0.034927, 0.42685, 0.42685];
P_h = 102.1445 * 133.32237;  % Pa
Lambda_M = 1.65 * 0.9;
Lambda_0 = 0.65;
S_basal = 0.862055 * 1.75;  % N/m
alpha_h = 47.569;  % degree

%beta_1 = 0.5; %0.5-0.9
beta_1=0.5;
%beta_2 = 0.05; %0.01-0.2; 0.01, 0.05, 0.2
beta_2 = 0.1
mu_1 = 0.2; %0.03-0.1; 0.01, 0.05, 0.2, 

ce = 657.89;
cc1 = 3091.6;
cc2 = 18.850;
cm1 = 103.34;
cm2 = 11.526;
cm3 = 0.7936;
MWH = 1;  % 4.981617e-22;  % Molecular Weight is taken to be 300 kDa = 4.981617e-22 kg
MWH1 = 0;
a_max = 350;
km1 = 0;
km2 = 0;
t = 0;
sign = 0;
vv = 0;
pp = 0;

rho = 1050.0;  % kg/m^3
PI = 3.14159265;
Tol = 0.00000000001;

% Constants
zeta_c = 1;
K_act = 0.1;
Kg_sh = 0;

% Homeostatic values
ro_h = 0.0003282;  % outer radius
h_h = 0.000021175;
Gc_h = 1.07;
Gm_h = 1.25;
Ge_1 = 2.19;
Ge_2 = 1.64;
stress_h = 120000.0;  % Pa
shear_h = 1.5;  % Pa

pk1_a = [];
pk2_a = [];
pk3_a = [];
pk1_p = [];
pk2_p = [];
pk3_p = [];
pm_a = [];
pm_p = [];
r_p = [];
pk_a = [];
T1_c = zeros(1, 3);
T2_c = zeros(1, 3);
T2_m = zeros(1, 3);
T1_c_old = 0;
T2_c_old = 0;
T2_m_old = 0;
error1 = 0;
error2 = 0;
alpha_R = zeros(1, 4);
alpha = zeros(1, 4);
alpha_a = 0;
L1 = 0;
L2 = 0;
L2_a = 0;
Lk2 = 0;
Lk2_a = 0;
Lk2_n = 0;
Le2_n1 = 0;
Le2_n2 = 0;
Lm2_n1 = 0;
Lm2_n2 = 0;
phi_c0 = 0;
dt = 0;
ri = 0;
ri_h = 0;
r_h = 0;
r = 0;
h = 0;
hm = 0;
hc = 0;
r_old = 0;
CI = zeros(1, 4);
CI_p = zeros(1, 4);
CIm = 0;
CIm_p = 0;
M = zeros(1, 4);
Mc = 0;
Mm = 0;
Me = 0;
total_M = 0;
Mc_h = 0;
Mm_h = 0;

ri_rih = 0;
h_hh = 0;
T_Th =0;
Sigma_Sigmah = 0;

a_mean = 0;
a_half = 0;
dummy1 = 0;
dummy2 = 0;
dummy3 = 0;
dummy4 = 0;
mP = zeros(1, 4);
mPm = 0;
mP_h = zeros(1, 4);
mPm_h = 0;
mP_p = zeros(1, 4);
mPm_p = 0;
shear = 0;
pressure = 0;
sigma_k = zeros(1, 4);
sigma_k_p = zeros(1, 4);
sigma_m = 0;
dwdL22 = 0;
ddwddL22 = 0;
dwdL22_m = 0;
dwdL22_c = 0;
ddwddL22_m = 0;
ddwddL22_c = 0;
dwdL12_c = 0;
fun_r = 0;
dfun_rdr = 0;
S = 0;
dSdL2 = 0;
t_act = 0;
dt_actdL2 = 0;
wt = 0;
time_step = 0;
final_time = 0;
num_pa = 0;
i = 0;
j = 0;
k = 0;
n = 0;
a = 0;
flag = 0;
num_it1 = 0;
num_it2 = 0;
check = 0;
Max_it = 2000;
string_1 = '';
string_2 = '';
fid1 = 0;
fid2 = 0;

k_2 = 0.001;
% mu_2d = 0.0;
Kg_mu = 0.01;
Kg = 2.0;
%Kg_sh = 52;
Kg_sh = 20;
km2 = Kg_sh / 3;

% Kg_mu = 11;
% km2 = 2;
% Kg = 1;
% Kg_sh = 50;
% km2 = 10;
% km2 = 40;
% Kg = 1.655;
% Kg_mu = 0.420;

fprintf('\nKg=%e, Kg_mu=%e\n', Kg, Kg_mu);

% Output file
% string_1 = sprintf('p4_Kg11%.3f_Kg_expo1S_mu%.3f_k2%.3f_mu_2%.3f_Kg_sh%i.txt', Kg, Kg_mu, k_2, mu_2, floor(Kg_sh));
% string_1 = sprintf('p4_k_2_%.3f_Kg%.3f_mu_2%0.3f.txt', k_2, Kg, mu_2);
string_1 = sprintf('Results_b1%.3f_b2%.3f_m1%.3f.txt', beta_1, beta_2, mu_1);
disp(string_1)
errmg='';
[fid1 errmg]= fopen(string_1, 'w');
if fid1 == -1
    disp(errmg);
    fprintf('Error in opening a file\n');
    return;
end

% string_2 = sprintf('p4_Kg11%.3f_Kg_mu%.3f_km%.3f_mu_2%.3f_Kg_sh%i.info', Kg, Kg_mu, k_2, mu_2, floor(Kg_sh));
% % string_2 = sprintf('p4_k_2_%.3f_Kg%.3f_mu_2%0.3f.info', k_2, Kg, mu_2);
string_2 = sprintf('Parameters_b1%.3f_b2%.3f_m1%.3f.txt', beta_1, beta_2, mu_1);
fid2 = fopen(string_2, 'w');
if fid2 == -1
    fprintf('Error in opening a file\n');
    return;
end

% Computational parameters
time_step = 20; % Number of time steps per day
final_time = 200; % Duration of time for computation (days)
% final_time = 56; % Duration of time for computation (days)
dt = 1.0 / time_step;
num_pa = time_step * a_max + 1;

% Axial direction is fixed
L1 = 1.0;

% Memory allocation
pk1_a = zeros(num_pa, 1);
pk2_a = zeros(num_pa, 1);
pk3_a = zeros(num_pa, 1);
pk1_p = zeros(num_pa, 1);
pk2_p = zeros(num_pa, 1);
pk3_p = zeros(num_pa, 1);
pm_a = zeros(num_pa, 1);
pm_p = zeros(num_pa, 1);
r_p = zeros(num_pa, 1);

% Calculation for homeostatic variables
total_M = h_h * rho;
phi_c0 = 1.0 - phi_e0 - phi_m0 - phi_f;
Mc_h = total_M * phi_c0;
Mm_h = total_M * phi_m0;
fprintf('Mc_h=%e\t total_M=%e\n', Mc_h, total_M);

ri_h = ro_h - h_h;
r_h = 0.5 * (ri_h + ro_h);

% Homeostatic values for pk, CI, a_mean, a_half
a_mean = 0;
pk1_p(1) = 1.0;
pk2_p(1) = 1.0;
pk3_p(1) = 1.0;
pm_p(1) = 1.0;
r_p(1) = r_h;

flag = 1;
t = 0.0;

for i = 2:num_pa
%     pk1_p(i) = pk1_p(i-1) * exp(-0.5*dt*(mu_F((i-1)*dt,0)+mu_F((i-2)*dt,0)));
%     pk2_p(i) = pk2_p(i-1) * exp(-0.5*dt*(mu_F((i-1)*dt,0)+mu_F((i-2)*dt,0)));
%     pk3_p(i) = pk3_p(i-1) * exp(-0.5*dt*(mu_F((i-1)*dt,0)+mu_F((i-2)*dt,0)));
%     pm_p(i) = pm_p(i-1) * exp(-0.5*dt*(mu_F((i-1)*dt,0)+mu_F((i-2)*dt,0)));  

    pk1_p(i) = pk1_p(i-1) * exp(-0.5*dt*(0.02));
    pk2_p(i) = pk2_p(i-1) * exp(-0.5*dt*(0.02));
    pk3_p(i) = pk3_p(i-1) * exp(-0.5*dt*(0.02));
    pm_p(i) = pm_p(i-1) * exp(-0.5*dt*(0.02));  
    
    % finding a half-life time
    if flag && (pk1_p(i) > 0.5)
        dummy1 = (i-1) * dt;
        dummy2 = i * dt;
        dummy3 = pk1_p(i-1);
        dummy4 = pk1_p(i);
        if dummy4 ~= dummy3
            a_half = dummy1 + dt * (0.5 - dummy3) / (dummy4 - dummy3);
        end
        flag = 0;
    end
    
    a_mean = a_mean + 0.5 * dt * (pk2_p(i-1) + pk2_p(i));
    r_p(i) = r_h;
end

mu_2= 1/a_mean;

fprintf('\nhalf-life time=%e\n', a_half);
fprintf('mean life time = %e\n', a_mean);

% beta_2d = k_2 / (k_2 + mu_2d);
% beta_2 = 0.05;

Mc = 0.0;
for i = 1:4
    M(i) = Mc_h * phi_k0(i);
    Mc = Mc + M(i);
    mP(i) = (beta_2+mu_1)*M(i)*mu_2/(MWH*beta_1*beta_2);
    mP_p(i) = mP(i);
    mP_h(i) = mP(i);
    CI(i) = M(i)*mu_2 / (MWH * beta_2);
    CI_p(i) = CI(i);
    sigma_k_p(i) = stress_h;
end

Me = total_M * phi_e0;
Mm = total_M * phi_m0;
CIm = Mm*mu_2 / (MWH * beta_2);
CIm_p = CIm;
mPm = (beta_2+mu_1)*Mm*mu_2/(MWH*beta_1*beta_2);
mPm_p = mPm;

% mPm_h = mPm_p;   % Commented out on October 14, 2015
mPm_h = Mm / a_mean;   % Added on October 14, 2015

alpha_R = [0.0, PI/2.0, alpha_h*PI/180.0, 2.0*PI - alpha_h*PI/180.0];

for i = 1:num_pa
    pk1_p(i) = pk1_p(i) * M(1) / a_mean;
    pk2_p(i) = pk2_p(i) * M(2) / a_mean;
    pk3_p(i) = pk3_p(i) * M(3) / a_mean;
    pm_p(i) = pm_p(i) * Mm / a_mean;
end

for i = 1:3
    T1_c(i) = stress_h;
    T2_c(i) = stress_h;
    T2_m(i) = stress_h;
end

t = 0.0;

% Save information of parameters
fprintf(fid2, 'The homeostatic value for stress in media includes active one\n');
fprintf(fid2, 'k_2=%e, mu_2=%e, beta_1=%e, zeta_c=%e, K_act=%e\n', k_2, mu_2, beta_1, zeta_c, K_act);
fprintf(fid2, 'Kg_mu=%e, Kg=%e, km2=%e, Kg_sh=%ee\n', Kg_mu, Kg, km2, Kg_sh);
fprintf(fid2, 'Homeostatics stress = %e\n', stress_h);
fprintf(fid2, 'Gch=%e, Gmh=%e, Gez=%e, Get=%e\n', Gc_h, Gm_h, Ge_1, Ge_2);
fprintf(fid2, 'L0=%e, Lm=%e, ro_h=%e, r_h=%e, h_h=%e\n', Lambda_0, Lambda_M, ro_h, r_h, h_h);
fprintf(fid2, 'ce=%e, cc1=%e, cc2=%e\n', ce, cc1, cc2);
fprintf(fid2, 'cm1=%e, cm2=%e, cm3=%e\n', cm1, cm2, cm3);
fprintf(fid2, 'phi_k0=[%e, %e, %e, %e]\n', phi_k0(1), phi_k0(2), phi_k0(3), phi_k0(4));
fprintf(fid2, 'P_h=%e, S_basal=%e, alpha_h=%e\n', P_h, S_basal, alpha_h);




% begin of a for-loop for the time
sum = 0;
r_act = r_h;
while (t <= final_time)
    t = t + dt;
    pressure = Given_P(t); % pascal
    dummy1 = abs((Given_P(t) - Given_P(t - dt)) / Given_P(t));
    dummy2 = abs((Given_flow(t) - Given_flow(t - dt)) / Given_flow(t));
    dummy3 = abs((Given_P(t - dt) - Given_P(t - 2.0 * dt)) / Given_P(t - dt));
    dummy4 = abs((Given_flow(t - dt) - Given_flow(t - 2.0 * dt)) / Given_flow(t - dt));
    
    % if the change in the flow is significant then predictor is set to the 
    % previous value, otherwise use 1st order predictor
    if dummy1 + dummy2 + dummy3 + dummy4 > 0.01
        r = r_p(1);
        T1_c(1) = T1_c(2);
        T2_c(1) = T2_c(2);
        T2_m(1) = T2_m(2);
    else
        T1_c(1) = 0.5 * (3.0 * T1_c(2) - T1_c(3));
        T2_c(1) = 0.5 * (3.0 * T2_c(2) - T2_c(3));
        T2_m(1) = 0.5 * (3.0 * T2_m(2) - T2_m(3));
        r = 0.5 * (3.0 * r_p(1) - r_p(2));
    end

    
    L2 = r / r_h;
    sum = sum + K_act * L2 * dt * exp(K_act * t); % added for r_act, Lm_act
    r_act = r_act + K_act * (r - r_act) * dt;


alpha(1) = 0.0;
alpha(2) = pi/2.0;
alpha(3) = atan(L2/L1*tan(alpha_h*pi/180));
alpha(4) = 2.0*pi - alpha(3);
zeta0 = 0; zeta1 = 0; zeta2 = 0; zetam = 0;

for a = 1:num_pa
    if a == 1
        L2_a = r / r_h;
        wt = 0.5 * dt;
    else
        L2_a = r_p(a-1) / r_h;
        if a == num_pa
            wt = 0.5 * dt;
        else
            wt = dt;
        end
    end
    
    for j = 1:3
        if j == 1 || j == 2
            alpha_a = alpha(j);
        end
        
        if j == 3
            alpha_a = atan(L2_a / L1 * tan(alpha_h * pi / 180));
        end
        
        Lk2_a = (L1 * cos(alpha_a))^2 + (L2_a * sin(alpha_a))^2;
        Lk2 = (L1 * cos(alpha(j)))^2 + (L2 * sin(alpha(j)))^2;
        
        if Lk2_a > Tol
            Lk2_n = Gc_h^2 * Lk2 / Lk2_a;
            if Lk2_n < 1
                Lk2_n = 1;
            end
            
            if j == 1
                zeta0 = zeta0 + pk1_a(a) * sqrt(Lk2_n) * (Lk2_n - 1) * exp(cc2 * (Lk2_n - 1)^2) * wt;
            elseif j == 2
                zeta1 = zeta1 + pk2_a(a) * sqrt(Lk2_n) * (Lk2_n - 1) * exp(cc2 * (Lk2_n - 1)^2) * wt;
            elseif j == 3
                zeta2 = zeta2 + pk3_a(a) * sqrt(Lk2_n) * (Lk2_n - 1) * exp(cc2 * (Lk2_n - 1)^2) * wt;
            end

        end
    end
    
    Lm2_n1 = Gm_h^2;
    
    if L2_a > Tol
        Lm2_n2 = Gm_h^2 * L2^2 / (L2_a^2);
        if Lm2_n2 < 1
            Lm2_n2 = 1;
        end
        
        zetam = zetam + pm_a(a) * Gm_h * (L2 / L2_a) * (cm1 * (1 - 1 / (Lm2_n1 * Lm2_n2^2)) + cm2 * (Lm2_n2 - 1) * exp(cm3 * (Lm2_n2 - 1)^2)) * wt;
    end
end


zeta0 = zeta0 / (M(1) * Gc_h * (Gc_h^2 - 1) * exp(cc2 * (Gc_h^2 - 1)^2));
zeta1 = zeta1 / (M(2) * Gc_h * (Gc_h^2 - 1) * exp(cc2 * (Gc_h^2 - 1)^2));
zeta2 = zeta2 / (M(3) * Gc_h * (Gc_h^2 - 1) * exp(cc2 * (Gc_h^2 - 1)^2));

zetam = 1;    % added on March 31, 2015, to make smooth muscle cells degradation rate constant

for i = 1:4
    sigma_k(i) = sqrt((T1_c(1) * cos(alpha(i)))^2 + (T2_c(1) * sin(alpha(i)))^2);
    
    if t <= dt
        shear = shear_h;
    end
    
    mP(i) = f_mP(Mc_h, Mc, sigma_k(i), mP_h(i), Kg, shear);
    CI(i) = CI_p(i) +   dt * (beta_1 * mP(i)- (beta_2+mu_1)*CI_p(i));
end

if MWH1 == 0
    MWH1 = MWH;
else
    MWH1 = L2 * L1 * h * MWH / h_h;
end

pk1_a = MWH1 * k_2 * CI(1);
pk2_a = MWH1 * k_2 * CI(2);
pk3_a = MWH1 * k_2 * CI(3);

mPm = f_mP(Mm_h, Mm, T2_m(1), mPm_h, Kg, shear);
CIm = CIm_p + dt * (beta_1 * mPm- (beta_2+mu_1)*CIm_p);
pm_a = mPm;

M(1) = 0.0;
M(2) = 0.0;
M(3) = 0.0;
Mm = 0.0;

% prediction for the age distribution function
for i = 2:num_pa
    pk1_a(i) = pk1_p(i-1) * exp(-0.5 * dt * (mu_F(i*dt, zeta0) + mu_F((i-1)*dt, zeta0)));
    pk2_a(i) = pk2_p(i-1) * exp(-0.5 * dt * (mu_F(i*dt, zeta1) + mu_F((i-1)*dt, zeta1)));
    pk3_a(i) = pk3_p(i-1) * exp(-0.5 * dt * (mu_F(i*dt, zeta2) + mu_F((i-1)*dt, zeta2)));
    pm_a(i) = pm_p(i-1) * exp(-0.5 * dt * (mu_F(i*dt, zetam) + mu_F((i-1)*dt, zetam)));
    
    M(1) = M(1) + 0.5 * dt * (pk1_a(i) + pk1_a(i-1));
    M(2) = M(2) + 0.5 * dt * (pk2_a(i) + pk2_a(i-1));
    M(3) = M(3) + 0.5 * dt * (pk3_a(i) + pk3_a(i-1));
    Mm = Mm + 0.5 * dt * (pm_a(i) + pm_a(i-1));
end

M(4) = M(3);
Mc = M(1) + M(2) + M(3) + M(4);
total_M = (Mc + Mm + Me) / (1 - phi_f);
h = total_M / (rho * L1 * L2);

% begin of while-loop for mass mp
num_it1 = 0;
num_it2 = 0;
error1 = 10.0;

while (error1 > Tol && (num_it1 < Max_it && num_it2 < Max_it))
    num_it1 = num_it1 + 1;
    
    if (num_it1 == Max_it)
        fprintf('Iteration reached maximum number !!\n');
        final_time = 0;
    end
 T1_c_old = T1_c(1);
T2_c_old = T2_c(1);
T2_m_old = T2_m(1);

% begin of while-loop for r(t) using a N-R method
num_it2 = 0;
error2 = 10.0;

while (error2 > Tol && num_it2 < Max_it)
    num_it2 = num_it2 + 1;
    
    if (num_it2 == Max_it)
        fprintf('Iteration for N-R method reached maximum number!!\n');
        final_time = 0;
    end
    
    r_old = r;
    L2 = r / r_h;
    h = total_M / (rho * L1 * L2);
    
    ri = r - 0.5 * h;
    shear = shear_h * Given_flow(t) * (ri_h / ri) ^ 3;
    
    alpha = zeros(1, 4);
    alpha(1) = 0.0;
    alpha(2) = pi / 2.0;
    alpha(3) = atan(L2 / L1 * tan(alpha_h * pi / 180));
    alpha(4) = 2.0 * pi - alpha(3);
    
    dwdL22_c = 0.0;
    dwdL12_c = 0.0;
    dwdL22_m = 0.0;
    ddwddL22_c = 0.0;
    ddwddL22_m = 0.0;
    
    Le2_n1 = Ge_1 * Ge_1 * L1 * L1;
    Le2_n2 = Ge_2 * Ge_2 * L2 * L2;
    dwdL22 = Me * (ce / 2.0) * Ge_2 * Ge_2 * (1.0 - 1.0 / (Le2_n1 * Le2_n2 * Le2_n2));
    ddwddL22 = Me * ce * (Ge_2 ^ 4) / (Le2_n1 * (Le2_n2 ^ 3));
    
    for a = 0:num_pa-1
    if a == 0
        L2_a = r / r_h;
        wt = 0.5 * dt;
    else
        L2_a = r_p(a) / r_h;
        if a == num_pa - 1
            wt = 0.5 * dt;
        else
            wt = dt;
        end
    end
    
    % Four fiber families
    for i = 0:3
        if i == 0 || i == 1
            alpha_a = alpha(i + 1);
        end
        if i == 2
            alpha_a = atan(L2_a / L1 * tan(alpha_h * pi / 180));
        end
        if i == 3
            alpha_a = 2.0 * pi - atan(L2_a / L1 * tan(alpha_h * pi / 180));
        end
        
        Lk2_a = L1 * cos(alpha_a) ^ 2 + L2_a * sin(alpha_a) ^ 2;
        Lk2 = L1 * cos(alpha(i + 1)) ^ 2 + L2 * sin(alpha(i + 1)) ^ 2;
        
        if Lk2_a > Tol
            Lk2_n = Gc_h * Gc_h * Lk2 / Lk2_a;
        end
        if Lk2_n < 1
            Lk2_n = 1;
        end
        
        if i == 0
            pk_a = pk1_a(a + 1);
        end
        if i == 1
            pk_a = pk2_a(a + 1);
        end
        if i == 2 || i == 3
            pk_a = pk3_a(a + 1);
        end
        
        dwdL22_c = dwdL22_c + 0.5 * cc1 * pk_a * Gc_h * Gc_h / Lk2_a * sin(alpha_a) ^ 2 ...
            * (Lk2_n - 1) * exp(cc2 * (Lk2_n - 1) * (Lk2_n - 1)) * wt;
        dwdL12_c = dwdL12_c + 0.5 * cc1 * pk_a * Gc_h * Gc_h / Lk2_a * cos(alpha_a) ^ 2 ...
            * (Lk2_n - 1) * exp(cc2 * (Lk2_n - 1) * (Lk2_n - 1)) * wt;
        ddwddL22_c = ddwddL22_c + 0.5 * cc1 * pk_a * Gc_h * Gc_h * Gc_h * Gc_h / (Lk2_a * Lk2_a) ...
            * sin(alpha_a) ^ 4 * (1.0 + 2.0 * cc2 * (Lk2_n - 1) * (Lk2_n - 1)) ...
            * exp(cc2 * (Lk2_n - 1) * (Lk2_n - 1)) * wt;
    end
    
    Lm2_n1 = Gm_h * Gm_h;
    Lm2_n2 = Gm_h * Gm_h * L2 * L2 / (L2_a * L2_a);
    
   if Lm2_n2 < 1
    Lm2_n2 = 1;
   end
   dwdL22_m = dwdL22_m + pm_a(a + 1) * Gm_h * Gm_h / (L2_a * L2_a) ...
    * (0.5 * cm1 * (1 - 1.0 / (Lm2_n1 * Lm2_n2 * Lm2_n2)) ...
    + 0.5 * cm2 * (Lm2_n2 - 1) * exp(cm3 * (Lm2_n2 - 1) * (Lm2_n2 - 1))) * wt;
   ddwddL22_m = ddwddL22_m + pm_a(a + 1) * (Gm_h / L2_a)^4 ...
    * (cm1 / (Lm2_n1 * Lm2_n2^3) ...
    + 0.5 * cm2 * (1 + 2 * cm3 * (Lm2_n2 - 1) * (Lm2_n2 - 1)) ...
    * exp(cm3 * (Lm2_n2 - 1) * (Lm2_n2 - 1))) * wt;
   
   end

 % Active muscle tone
 h = total_M / rho;
hm = Mm / ((1 - phi_f) * rho * L1 * L2);
sigma_m = 2 * L2 * dwdL22_m / (L1 * hm);

S = (1 / (1 - exp(-km2^2 / 400))) * (Mm / Mm_h) * S_basal * (1 - exp(-((km2 / 20 - km2 * (shear / shear_h - 1))^2)));

L2m_act = r / r_act;


dsdr = 6 * (1 / (1 - exp(-km2^2 / 400))) * (Mm / Mm_h) * S_basal ...
    * (km2 / 20 - km2 * (shear / shear_h - 1)) * (km2 * shear / (ri * shear_h)) ...
    * (1 + 0.5 * total_M / (rho * L1 * L2 * L2 * r_h)) ...
    * exp(-((km2 / 20 - km2 * (shear / shear_h - 1))^2));

% if S > 4 * S_basal
%     S = 4 * S_basal;
%     dsdr = 0;
% end

if ((Lambda_M - L2m_act) / (Lambda_M - Lambda_0))^2 <= 1
    t_act = S * L2m_act * (1 - ((Lambda_M - L2m_act) / (Lambda_M - Lambda_0))^2);

    if r ~= 0
        dt_actdr = dsdr * L2m_act * (1 - ((Lambda_M - L2m_act) / (Lambda_M - Lambda_0))^2) + t_act / r + 2 * L2m_act * S * (Lambda_M - L2m_act) * L2m_act / (r * (Lambda_M - Lambda_0)^2);
    end
else
    t_act = 0;
    dt_actdr = 0;
end

dwdL22 = dwdL22 + dwdL22_m + dwdL22_c;
ddwddL22 = ddwddL22 + ddwddL22_m + ddwddL22_c;
hc = Mc / ((1 - phi_f) * rho * L1 * L2);
fun_r = 2.0 * L2 * dwdL22 / L1 + t_act - pressure * r;
dfun_rdr = (1.0 / r_h) * (2 * dwdL22 / L1 + 4 * L2 * L2 * ddwddL22 / L1) + dt_actdr - pressure;

if abs(dfun_rdr) < Tol
    if dfun_rdr < 0
        dfun_rdr = -Tol;
    else
        dfun_rdr = Tol;
    end
end

r = r - 0.2 * fun_r / dfun_rdr;

r_act = r_act + K_act * (r - r_old);

if r_old ~= 0
    error2 = abs((r_old - r) / r_old);
end

end

% End of while-loop for r(t)
hm = Mm / ((1 - phi_f) * rho * L1 * L2);
hc = Mc / ((1 - phi_f) * rho * L1 * L2);
T1_c(1) = 2 * L1 * dwdL12_c / (L2 * hc);
T2_c(1) = 2 * L2 * dwdL22_c / (L1 * hc);
T2_m(1) = (2 * L2 * dwdL22_m / L1 + t_act) / hm;
L2 = r / r_h;

% Commented on June 20, 2014.
alpha(2) = atan(L2 / L1 * tan(alpha_h * pi / 180));
alpha(3) = 2.0 * pi - alpha(2);

zeta0 = 0;
zeta1 = 0;
zeta2 = 0;
zetam = 0;

for a = 1:num_pa
    if a == 1
        L2_a = r / r_h;
        wt = 0.5 * dt;
    else
        L2_a = r_p(a - 1) / r_h;
        
        if a == num_pa
            wt = 0.5 * dt;
        else
            wt = dt;
        end
    end
    
    for j = 0:2
        if j == 0 || j == 1
            alpha_a = alpha(j + 1);
        end
        
        if j == 2
            alpha_a = atan(L2_a / L1 * tan(alpha_h * pi / 180));
        end
        
        Lk2_a = L1^2 * cos(alpha_a)^2 + L2_a^2 * sin(alpha_a)^2;
        Lk2 = L1^2 * cos(alpha(j + 1))^2 + L2^2 * sin(alpha(j + 1))^2;
        
        if Lk2_a > Tol
            Lk2_n = Gc_h^2 * Lk2 / Lk2_a;
        end
        
        if Lk2_n < 1
            Lk2_n = 1;
        end
        
        if j == 0
            zeta0 = zeta0 + pk1_a(a) * sqrt(Lk2_n) * (Lk2_n - 1) * exp(cc2 * (Lk2_n - 1)^2) * wt;
        elseif j == 1
            zeta1 = zeta1 + pk2_a(a) * sqrt(Lk2_n) * (Lk2_n - 1) * exp(cc2 * (Lk2_n - 1)^2) * wt;
        elseif j == 2
            zeta2 = zeta2 + pk3_a(a) * sqrt(Lk2_n) * (Lk2_n - 1) * exp(cc2 * (Lk2_n - 1)^2) * wt;
        end
    end
    
    Lm2_n1 = Gm_h^2;
    
    if L2_a > Tol
        Lm2_n2 = Gm_h^2 * L2^2 / (L2_a^2);
        
        if Lm2_n2 < 1
            Lm2_n2 = 1;
        end
        
        zetam = zetam + pm_a(a) * Gm_h * (L2 / L2_a) * (cm1 * (1 - 1 / (Lm2_n1 * Lm2_n2^2)) + cm2 * (Lm2_n2 - 1) * exp(cm3 * (Lm2_n2 - 1)^2)) * wt;
    end
end


zeta0 = zeta0 / (M(1) * Gc_h * (Gc_h^2 - 1) * exp(cc2 * (Gc_h^2 - 1)^2));
zeta1 = zeta1 / (M(2) * Gc_h * (Gc_h^2 - 1) * exp(cc2 * (Gc_h^2 - 1)^2));
zeta2 = zeta2 / (M(3) * Gc_h * (Gc_h^2 - 1) * exp(cc2 * (Gc_h^2 - 1)^2));
zetam = 1; % added on March 31, 2015, to make smooth muscle cells degradation rate constant

for i = 0:3
    sigma_k(i + 1) = sqrt(T1_c(1) * cos(alpha(i + 1))^2 + T2_c(1) * sin(alpha(i + 1))^2);
    mP(i + 1) = f_mP(Mc_h, Mc, sigma_k(i + 1), mP_h(i + 1), Kg, shear);

    CI(i + 1) = (CI_p(i + 1) + dt * (beta_1 * mP_p(i + 1) - (beta_2 + mu_1) * CI_p(i + 1))) / (1.0 + dt * (beta_2+mu_1));
end
pk1_a = MWH1 * k_2 * CI(1);
pk2_a = MWH1 * k_2 * CI(2);
pk3_a = MWH1 * k_2 * CI(3);
mPm = f_mP(Mm_h, Mm, T2_m(1), mPm_h, Kg, shear);
CIm = (CIm_p + dt * (beta_1 * mPm_p - (beta_2 + mu_1) * CIm_p)) / (1.0 + dt * (beta_2+mu_1));
pm_a = mPm;
M(1) = 0;
M(2) = 0;
M(3) = 0;
Mm = 0;

for i = 2:num_pa
    pk1_a(i) = pk1_p(i - 1) * exp(-0.5 * dt * (mu_F(i * dt, zeta0) + mu_F((i - 1) * dt, zeta0)));
    pk2_a(i) = pk2_p(i - 1) * exp(-0.5 * dt * (mu_F(i * dt, zeta1) + mu_F((i - 1) * dt, zeta1)));
    pk3_a(i) = pk3_p(i - 1) * exp(-0.5 * dt * (mu_F(i * dt, zeta2) + mu_F((i - 1) * dt, zeta2)));
    pm_a(i) = pm_p(i - 1) * exp(-0.5 * dt * (mu_F(i * dt, zetam) + mu_F((i - 1) * dt, zetam)));
    M(1) = M(1) + 0.5 * dt * (pk1_a(i) + pk1_a(i - 1));
M(2) = M(2) + 0.5 * dt * (pk2_a(i) + pk2_a(i - 1));
M(3) = M(3) + 0.5 * dt * (pk3_a(i) + pk3_a(i - 1));
Mm = Mm + 0.5 * dt * (pm_a(i) + pm_a(i - 1));
end

M(4) = M(3);
Mc = M(1) + M(2) + M(3) + M(4);
total_M = (Mc + Mm + Me) / (1 - phi_f);
h = total_M / (rho * L1 * L2);

if (T1_c_old ~= 0 && T2_c_old ~= 0 && T2_m_old ~= 0)
    error1 = abs((T1_c_old - T1_c(1)) / T1_c_old) + abs((T2_c_old - T2_c(1)) / T2_c_old) + abs((T2_m_old - T2_m(1)) / T2_m_old);
end
end 



for i = 1:4
    mP_p(i) = mP(i);
    CI_p(i) = CI(i);
    sigma_k_p(i) = sigma_k(i);
end
mPm_p = mPm;
CIm_p = CIm;

for i = 1:num_pa
    pk1_p(i) = pk1_a(i);
    pk2_p(i) = pk2_a(i);
    pk3_p(i) = pk3_a(i);
    pm_p(i) = pm_a(i);
    
    if i == num_pa
        r_p(1) = r;
    else
        r_p(num_pa - i + 1) = r_p(num_pa - i);
    end
end

for i = 1:2
    T1_c(i) = T1_c(i + 1);
    T2_c(i) = T2_c(i + 1);
    T2_m(i) = T2_m(i + 1);
end

T1_c(3) = 0.5 * (3 * T1_c(2) - T1_c(1));
T2_c(3) = 0.5 * (3 * T2_c(2) - T2_c(1));
T2_m(3) = 0.5 * (3 * T2_m(2) - T2_m(1));

ri_rih = ri/ri_h;
h_hh = h/h_h;

% re-setting shear_h and stress_h at the end of homeostatic condition period

if (t>1.5 && t<1.9)
    t
    shear_h = shear
    stress_h = T2_c(2)
end

T_Th =shear/shear_h;
Sigma_Sigmah = T2_c(2)/stress_h;

fprintf('time=%.2f, r=%e \t%e\n', t, r, zeta1);
fprintf(fid1, '%4.3f\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n', t, r, h, T1_c(1), T2_c(1), T2_m(1), shear, pressure * r / h, t_act / hm, S, L2m_act, Mm, ri_rih, h_hh, T_Th, Sigma_Sigmah);




   end

fclose(fid1);
fclose(fid2);
clear pk1_a pk2_a pk3_a pk1_p pk2_p pk3_p pm_a pm_p r_p;
disp('end');




 % Intramural Pressure as a function of time (day)
   function P = Given_P(time)
    P_h = 102.1445 * 133.32237;  % Pa
    if time >= 2
        P = 1.0 * P_h;
    else
        P = P_h;
    end
end

% Relative increase in flow
function flow = Given_flow(time)
    if time >= 2
        flow = 1.1;
    else
        flow = 1.0;
    end
end

% Rate of degradation
function rate = mu_F(age, zeta)
   global zeta_c;
    Kg_mu = 0;  % Need to assign a value for Kg_mu
    if zeta < zeta_c
        rate = 1 / 100.0;
    else
        rate = 1 / 100 + power(zeta - zeta_c, 2) * Kg_mu / 100;
    end
end

% Calculate f_mP
function mP = f_mP(Mi_0, Mi_t, S, mP_b, kkgg, tau)
   global stress_h;
   global Kg_sh;
   global shear_h;
    mP = mP_b * (1.0 + kkgg * (S / stress_h - 1.0) - Kg_sh * (tau / shear_h - 1));
    if mP > 0
        mP = mP;
    else
        mP = 0.0;
    end
end

% Calculate absolute value
function absVal = abss(ab)
    if ab < 0
        absVal = -ab;
    else
        absVal = ab;
    end
end



