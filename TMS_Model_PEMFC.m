% =========================================================================
% Mission Thermal Management System Simulation
% NOTE:Comments explain intent, units, and data flow for faster review.
% =========================================================================

clear; clc;
A_FC_eff=23880542.797;%295m2   % Effective FC active area [m^2]
n_comp = 4;                    % Number of serial components/HXs: [HTS, DC/DC, DC/AC, Motor]
nx=5;                          % Axial nodes per HX (method of lines)
t_step=30;                     % Samples per mission step for time-expanded storage/plots

%% 1. Load data from Excel
% Expected columns: time_sec_, altitude_m_, speed_ms_, P_W_
data = readtable('Mission_Profile_Data_Total.xlsx');

%% 2. Extract parameters
time = data.time_sec_;         % [s]
%time_step = length(time);     % (unused)
altitude = data.altitude_m_;   % [m]
speed = data.speed_ms_;        % [m/s]
P_prop = data.P_W_;            % [W] propulsion power demand
P_prop_initial=P_prop(1);      % for normalized plotting
n = length(time);              % number of mission time steps: 933

%% 3. Preallocate arrays for results
% IMPORTANT: Temperatures/T_Comp are time-expanded (t_step samples per mission step).
%            
Powers = zeros(n, 8);          % [W] columns: [P_M, P_DCAC, P_HTS, P_DCDC, P_FC_el, P_comp, P_TMS, P_H2_LHV_eff]
Qtransfer = zeros(n, 5);       % [W] columns: [Q_M, Q_DCAC, Q_HTS, Q_DCDC, Q_FC_H2]
Temperatures = zeros(t_step*n, n_comp); % H2 outlet temps per HX 
T_Comp=zeros(t_step*n,n_comp);          % He/component temps per HX 
Temp_He=zeros(n,n_comp);                  % Placeholder 
mass = zeros(n, 2);            % [kg/s] [m_dot_H2, m_dot_air]
j = zeros(n, 1);               % [A/cm^2] FC current density
RE_HTS=zeros(n,1);             % Reynolds number of H2 leaving HTS HX (diagnostic)
Tnew= zeros(t_step,n_comp);    % per-step H2 temp block (interpolated)
T2new= zeros(t_step,n_comp);   % per-step He/comp temp block (interpolated)

%% 4. Initial fuel cell power and hydrogen temperature
P_FC_el_initial = 1.8e7;       % [W] initial guess for FC electrical power
P_star_TMS = 5e6;              % [W] initial guess for TMS power

% First mission point
ho = altitude(1);
so = speed(1);
po = P_prop(1);

k=0; % rolling row index into the expanded temperature matrices
[P,Q_dot,T,m_dot,jcell,y,Re_HTS] = PowerModel2j_IZEA_2_Initial(po, ho, so, A_FC_eff, P_FC_el_initial,P_star_TMS);

% Unpack first-step power split/state
P_FC_el = P(5);                % [W] FC electric power for next step
P_TMS = P(7);                  % [W] TMS power for next step
Powers(1,:) = P;               % power channels
Qtransfer(1,:) = Q_dot;        % heat loads

% Extract He/component temps at HX outlets from state 'y' (last axial node per HX)
T2=y(:,[nx*(n_comp+1),nx*(n_comp+2),nx*(n_comp+3),nx*(n_comp+4)]);

% Interpolate outlet temps within this mission step to 't_step' samples
for z=1:n_comp
    Tfix=linspace(T(1,z),T(end,z),t_step);     % H2 outlet per HX
    T2fix=linspace(T2(1,z),T2(end,z),t_step);  % He/comp outlet per HX
    Tnew(:,z)=(Tfix)';
    T2new(:,z)=(T2fix)';
end

% Append into expanded matrices
[m, ~] = size(Tnew); % current block size (=t_step)
m=k+m;               % end index for this append
Temperatures(k+1:m,:)=Tnew;
T_Comp(k+1:m,:)=T2new;
k=m;                 % advance rolling index

% Store first-step scalars
mass(1,:) = m_dot;
j(1) = jcell;
RE_HTS(1)=Re_HTS;

%% 5. Mission Iteration - consider vectorizing if PowerModel2j_IZEA_2_Optimized supports it
for i = 2:(n-1)
    fprintf("i= %.2f \n", i );
    % Current mission point
    h = altitude(i);
    s = speed(i);
    p = P_prop(i);
    
    % Power/thermal solve (warm-started with previous 'y', P_FC_el, P_TMS)
    tic;
    [P,Q_dot,T,m_dot,jcell,y,Re_HTS] = PowerModel2j_IZEA_2_Optimized(p, h, s, A_FC_eff, P_FC_el,P_TMS,y);
    toc;
    
    % He/comp outlet temps
    T2=y(:,[nx*(n_comp+1),nx*(n_comp+2),nx*(n_comp+3),nx*(n_comp+4)]);

    % Interpolate per-step outlets to 't_step' samples for smooth plots
    for z=1:n_comp
        Tfix=linspace(T(1,z),T(end,z),t_step);
        T2fix=linspace(T2(1,z),T2(end,z),t_step);
        Tnew(:,z)=(Tfix)';
        T2new(:,z)=(T2fix)';
    end

    % Store mission-time series
    Powers(i,:) = P;
    Qtransfer(i,:) = Q_dot;

    % Append expanded temps
    [m, ~] = size(Tnew);
    m=k+m;
    Temperatures(k+1:m,:)=Tnew;
    T_Comp(k+1:m,:)=T2new;
    k=m;

    % Scalars
    mass(i,:) = m_dot;
    j(i) = jcell;
    RE_HTS(i)=Re_HTS;

    % Update guesses for next step
    P_FC_el = P(5);
    P_TMS = P(7);
end

%% 6. Optional: Handle the last iteration if needed
% Replicate last valid row for mission-time arrays so length == n
if n > 1
    Powers(n,:) = Powers(n-1,:);
    Qtransfer(n,:) = Qtransfer(n-1,:);
    mass(n,:) = mass(n-1,:);
    j(n) = j(n-1);
    RE_HTS(n)=RE_HTS(n-1);
    Temp_He(n,:)=Temp_He(n-1,:);
end

%% Mission Profile Plot (normalized power, altitude, speed)
figure(1)
subplot(3,1,1)
plot(time/3600,P_prop/P_prop_initial,LineWidth=2,Color=[1 0 0]);
hold on
set(gca, 'LineWidth', 2);
set(gca,"FontSize",12,"FontWeight","bold")
yticks([0 0.2 0.4 0.6 0.8 1]);
ylabel("P[u]","FontSize",15);
grid on
grid(gca, 'on');
set(gca, 'GridAlpha', 0.5); % styling
set(gca, 'GridLineStyle', '-');
set(gca, 'GridColor', 'k');
set(gca, 'GridLineWidth', 1);
hold off

subplot(3,1,2)
plot(time/3600,altitude/1000,LineWidth=2,Color=[0 1 0]);
hold on
set(gca, 'LineWidth', 2);
set(gca,"FontSize",12,"FontWeight","bold")
yticks([0 4 8 12 16]); % manual ticks for readability
ylabel("h[km]","FontSize",15);
grid on
grid(gca, 'on');
set(gca, 'GridAlpha', 0.5);
set(gca, 'GridLineStyle', '-');
set(gca, 'GridColor', 'k');
set(gca, 'GridLineWidth', 1);
hold off

subplot(3,1,3)
plot(time/3600,speed*3.6,LineWidth=2,Color=[0 0 1]);
hold on
set(gca, 'LineWidth', 2);
set(gca,"FontSize",12,"FontWeight","bold")
yticks([0 200 400 600 800]);
xlabel("time [h]","FontSize",15);
ylabel("v[km/h]","FontSize",15);
grid on
grid(gca, 'on');
set(gca, 'GridAlpha', 0.5);
set(gca, 'GridLineStyle', '-');
set(gca, 'GridColor', 'k');
set(gca, 'GridLineWidth', 1);
hold off

%% save the figure plot
print('mission_profile.jpg','-djpeg');

%% Systems Efficiency (proxy)
% Using P_H2_LHV_eff = P_FC_el / eff_FC(load) as FC input proxy → P_prop / P_H2_LHV_eff
P_H2_LHV_eff = Powers(:,8);
efficiency = P_prop./P_H2_LHV_eff;

%% Plot Important Parameters (2x2 grid)
figure(4)
subplot(2,2,1)
plot(time/3600,P_prop/1000000,LineWidth=2,Color=[1 0 0]);
hold on
set(gca, 'LineWidth', 2);
title("Power Prop","FontSize",20,"FontWeight","bold");
set(gca,"FontSize",15,"FontWeight","bold")
xticks([2 4 6 8]);
ylabel("P[MW]","FontSize",15,"FontWeight","bold");
grid on
grid(gca, 'on');
set(gca, 'GridAlpha', 0.5);
set(gca, 'GridLineStyle', '-');
set(gca, 'GridColor', 'k');
set(gca, 'GridLineWidth', 1);
hold off

subplot(2,2,2)
plot(time/3600,j,LineWidth=2,Color=[0 1 0]);   
hold on
set(gca, 'LineWidth', 2);
title("Current density","FontSize",20,"FontWeight","bold");
xticks([2 4 6 8]);
ylabel("j[A/cm2]","FontSize",15,"FontWeight","bold");
set(gca,"FontSize",15,"FontWeight","bold")
grid on
grid(gca, 'on');
set(gca, 'GridAlpha', 0.5);
set(gca, 'GridLineStyle', '-');
set(gca, 'GridColor', 'k');
set(gca, 'GridLineWidth', 1);
hold off

subplot(2,2,3)
plot(time/3600,Powers(:,6)/1000000,LineWidth=2,Color=[0 0 1]);
hold on
set(gca, 'LineWidth', 2);
title("Compressor Power","FontSize",20,"FontWeight","bold");
xlabel("Time [h]","FontSize",15,"FontWeight","bold");
ylabel("P[MW]","FontSize",15,"FontWeight","bold");
xticks([2 4 6 8]);
set(gca,"FontSize",15,"FontWeight","bold")
grid on
grid(gca, 'on');
set(gca, 'GridAlpha', 0.5);
set(gca, 'GridLineStyle', '-');
set(gca, 'GridColor', 'k');
set(gca, 'GridLineWidth', 1);
hold off

subplot(2,2,4)
plot(time/3600,efficiency*100,LineWidth=2,Color=[0 0 0]);
hold on
set(gca, 'LineWidth', 2);
title("Effciency","FontSize",20,"FontWeight","bold"); % (typo in label kept as-is)
xlabel("Time [h]","FontSize",15,"FontWeight","bold");
ylabel("η[%]","FontSize",15,"FontWeight","bold");
xticks([2 4 6 8]);
set(gca,"FontSize",15,"FontWeight","bold")
grid on
grid(gca, 'on');
set(gca, 'GridAlpha', 0.5);
set(gca, 'GridLineStyle', '-');
set(gca, 'GridColor', 'k');
set(gca, 'GridLineWidth', 1);
hold off

%% save the figure plot
print('important_parameters.jpg','-djpeg');
saveas(gcf, 'important_parameters.png'); % Duplicate save as PNG

%% P_FC profile (fuel cell electrical power)
figure(2);
plot(time/3600,Powers(:,5)/1000000,LineWidth=2,Color=[1 0 0]);
hold on
set(gca, 'LineWidth', 2);
title("Power FC","FontSize",20,"FontWeight","bold");
set(gca,"FontSize",15,"FontWeight","bold")
xlabel("Time [h]");
ylabel("Power [MW]","FontSize",15,"FontWeight","bold");
grid on
grid(gca, 'on');
set(gca, 'GridAlpha', 0.5);
set(gca, 'GridLineStyle', '-');
set(gca, 'GridColor', 'k');
set(gca, 'GridLineWidth', 1);
hold off

%% Reynolds Hydrogen leaving HTS HX (diagnostic)
figure(3);
plot(time/3600,RE_HTS,LineWidth=2,Color=[1 0 0]);
hold on
set(gca, 'LineWidth', 2);
title("Reynolds HTS","FontSize",20,"FontWeight","bold");
set(gca,"FontSize",15,"FontWeight","bold")
xlabel("Time [h]");
ylabel("RE","FontSize",15,"FontWeight","bold");
grid on
grid(gca, 'on');
set(gca, 'GridAlpha', 0.5);
set(gca, 'GridLineStyle', '-');
set(gca, 'GridColor', 'k');
set(gca, 'GridLineWidth', 1);
hold off
exportgraphics(gcf,'Reynolds_H2.jpg','Resolution',300);

%% Mass Flow (H2 and Air)
mass2 = mass';                  % rows: fluids, cols: time
fluid_names = {'H_2', 'Air'};  % labels for subplots
num_fluids = size(mass2, 1);

figure(7); clf
set(gcf, 'Position', [100, 100, 800, 1200]);  % Taller figure scaffolding

for i = 1:num_fluids
    subplot(num_fluids, 1, i)
    plot(time/3600, mass2(i, :), 'k', 'LineWidth', 2);
    grid on;
    ylabel('[kg/s]', 'FontSize', 12, 'FontWeight', 'bold');
    title(fluid_names{i}, 'FontSize', 14, 'FontWeight', 'bold');
    
    % Styling
    set(gca, 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'LineWidth', 1.2);
    set(gca, 'GridAlpha', 0.4);
    set(gca, 'GridColor', 'k');
    set(gca, 'GridLineStyle', '-');
    set(gca, 'GridLineWidth', 0.75);
    
    if i < num_fluids
        set(gca, 'XTickLabel', []);  % hide x labels except bottom
    else
        xlabel("Time [h]", 'FontSize', 14, 'FontWeight', 'bold');
    end
end

% Shared title and save
sgtitle('Mass Flow of Fluids Over Flight Mission', 'FontSize', 16, 'FontWeight', 'bold');
print('massflow.jpg','-djpeg','-r300');

%% Temperatures of hydrogen (H2 outlet per HX, expanded timeline)
% Matrix columns must match the component order used in TempCalc*: [HTS, DC/DC, DC/AC, Motor]
temp = Temperatures(:, 1:n_comp);
[r, ~] = size(temp);
t = linspace(0, n * t_step, r);        % synthetic time base [s] = (#steps * 30)
temp(end, :) = temp(end - 1, :);       % avoid trailing zero/NaN
temp2 = temp';

% Palettes and labels (order matches temp columns)
color_palette = lines(n_comp);
labels = {'HTS cable','DC/DC converter', 'DC/AC inverter', 'M motor'};

figure(5); clf
hold on
for i = 1:n_comp
    plot(t / 3600, temp2(i, :), 'LineWidth', 2, 'Color', color_palette(i, :));
end
xlabel("Time [h]", "FontSize", 20, "FontWeight", "bold");
ylabel("T [K]", "FontSize", 20, "FontWeight", "bold");
set(gca, "FontSize", 15, "FontWeight", "bold", "LineWidth", 1.2);
grid on;
set(gca, 'GridAlpha', 0.4, 'GridLineStyle', '-', 'GridColor', 'k', 'GridLineWidth', 0.75);
legend(labels, 'Location', 'eastoutside');
title('Temperatures of H2 at each HX outlet over time', 'FontSize', 18, 'FontWeight', 'bold');
hold off
print('temperatures.jpg','-djpeg','-r300');

%% Temperatures Helium after absorbing heat load (He/component outlets, expanded)
tempC = T_Comp(:, 1:n_comp);
[r, ~] = size(tempC);
t2 = linspace(0, n * t_step, r);       % synthetic time base [s]
tempC(end, :) = tempC(end - 1, :);     % avoid trailing zero/NaN
tempC2 = tempC';

figure(6); clf
hold on
for i = 1:n_comp
    plot(t2 / 3600, tempC2(i, :), 'LineWidth', 2, 'Color', color_palette(i, :));
end
xlabel("Time [h]", "FontSize", 20, "FontWeight", "bold");
ylabel("T [K]", "FontSize", 20, "FontWeight", "bold");
set(gca, "FontSize", 15, "FontWeight", "bold", "LineWidth", 1.2);
grid on;
set(gca, 'GridAlpha', 0.4, 'GridLineStyle', '-', 'GridColor', 'k', 'GridLineWidth', 0.75);
legend(labels, 'Location', 'eastoutside');
title('Temperature of Helium/Component over time', 'FontSize', 18, 'FontWeight', 'bold');
hold off
print('temperatures2.jpg','-djpeg','-r300');

% ============================== FUNCTIONS ===============================

function [P,Q_dot,T,m_dot,jcell,ys,Re_HTS] = PowerModel2j_IZEA_2_Optimized(Pprop,altitude,speed,A_FC_eff,P_FC_el_initial,P_star_TMS,y)
% Solve for P_FC_el and P_TMS (fsolve) at a mission point using warm-start,
% then compute component powers, heat loads, temps, mass flows, and Reynolds.
% Outputs:
%   P      [1x8]  = [P_M, P_DCAC, P_HTS, P_DCDC, P_FC_el, P_comp, P_TMS, P_H2_LHV_eff]
%   Q_dot  [1x5]  = [Q_M, Q_DCAC, Q_HTS, Q_DCDC, Q_FC_H2]
%   T      [t x 4]= H2 outlet temps [HTS, DC/DC, DC/AC, Motor]
%   m_dot  [1x2]  = [m_dot_H2, m_dot_air]
%   jcell  [A/cm^2]
%   ys           = ODE state history (for next warm-start)
%   Re_HTS       = Reynolds of H2 at HTS outlet

    % Load-percentage for efficiency fits (relative to 16.242 MW baseline)
    Pprop_initial = 16242000;
    load_percent_ratio = Pprop / Pprop_initial;
    
    % Quadratic fits for efficiencies vs load fraction
    eff_FC = [0.484 0.55 0.572 0.646];
    load_percent = [1, 0.7, 0.55, 0.30];
    eff_fc_fit = polyfit(load_percent, eff_FC, 2);
    eff_DCDC = [0.98 0.981 0.982 0.982];
    eff_DCDC_fit = polyfit(load_percent, eff_DCDC, 2);
    eff_HTS = [0.999 0.999 0.999 0.999];
    eff_HTS_fit = polyfit(load_percent, eff_HTS, 2);
    eff_motorD = [0.99 0.991 0.992 0.992];
    eff_motorD_fit = polyfit(load_percent, eff_motorD, 2);
    eff_motor = [0.992 0.994 0.995 0.997];
    eff_motor_fit = polyfit(load_percent, eff_motor, 2);
    
    % Evaluate fits
    v.eff_FC = polyval(eff_fc_fit, load_percent_ratio);
    v.eff_DCDC = polyval(eff_DCDC_fit, load_percent_ratio);
    v.eff_HTS = polyval(eff_HTS_fit, load_percent_ratio);
    v.eff_DCAC = polyval(eff_motorD_fit, load_percent_ratio);
    v.eff_M = polyval(eff_motor_fit, load_percent_ratio);
    
    % Downstream chain to get losses per block
    components = simPar_IZEA_2_Optimized(v, Pprop);
    P_M=components(1);
    Q_dot_M=components(2);
    P_DCAC=components(3);
    Q_dot_DCAC=components(4);
    P_HTS=components(5);
    Q_dot_HTS=components(6);
   
    % Algebraic system for [P_FC_el; P_TMS]
    target_fun = @(x) [
        x(1) - ((P_HTS + x(2) + PEMFCj(x(1), altitude, speed, A_FC_eff)) * ((1 / v.eff_DCDC) - 1)) ...
              - (P_HTS + x(2) + PEMFCj(x(1), altitude, speed, A_FC_eff));
        x(2) - 0.04 * 0.95 * ((x(1) / v.eff_FC) - x(1)) % 4% of FC heat; 95% to TMS
    ];

    % Initial guess
    x0 = [P_FC_el_initial, P_star_TMS];

    % fsolve with robust options + retry on failure (keeps sim running)
    opts = optimoptions('fsolve', ...
        'Display','off', ...
        'Algorithm','trust-region-dogleg', ...
        'StepTolerance',1e-9, ...
        'FunctionTolerance',1e-9, ...
        'MaxIterations',400, ...
        'MaxFunctionEvaluations',4000);

    [x_solution, fval, exitflag, output] = fsolve(target_fun, x0, opts);

    if exitflag <= 0 || any(~isfinite(x_solution))
        warning(['PowerModel2j_IZEA_2_Optimized: fsolve did not converge. ', ...
                 'exitflag=%d, iters=%d, normF=%g. Retrying with Levenberg–Marquardt.'], ...
                 exitflag, output.iterations, norm(fval));

        optsLM = optimoptions(opts, 'Algorithm','levenberg-marquardt');
        [x_solution2, fval2, exitflag2, output2] = fsolve(target_fun, x0, optsLM);

        if exitflag2 > 0 && all(isfinite(x_solution2))
            x_solution = x_solution2;  fval = fval2;  exitflag = exitflag2;  output = output2;
        else
            warning(['PowerModel2j_IZEA_2_Optimized: retry also failed. ', ...
                     'Keeping previous guess to continue simulation.']);
            x_solution = x0;  % safe fallback
        end
    end

    % Extract
    P_FC_el = x_solution(1);
    P_TMS = x_solution(2);

    % Compressor/DC-DC updates
    P_comp= PEMFCj_Optimized(P_FC_el, altitude, speed, A_FC_eff);
    Q_dot_DCDC = (P_HTS + P_comp+P_TMS) * ((1/v.eff_DCDC) - 1);
    P_DCDC = P_FC_el - Q_dot_DCDC;
    
    % Transient thermal (warm-started with 'y')
    [temperatures,m_dot,Q_dot_FC_H2, ys,Re_HTS]= TempCalcj_IZEA_2_Optimized(Pprop, P_FC_el, A_FC_eff, Q_dot_DCDC,v,y);
    
    % Pack outputs
    P_H2_LHV_eff = P_FC_el ./ v.eff_FC;
    P = [P_M, P_DCAC, P_HTS,P_DCDC, P_FC_el, P_comp,P_TMS, P_H2_LHV_eff];
    Q_dot = [Q_dot_M, Q_dot_DCAC, Q_dot_HTS,Q_dot_DCDC, Q_dot_FC_H2];
    T = temperatures;
    jcell = jcalculation(P_FC_el,A_FC_eff);
end

function [P,Q_dot,T,m_dot,jcell,yinitial,Re_HTS] = PowerModel2j_IZEA_2_Initial(Pprop,altitude,speed,A_FC_eff,P_FC_el_initial,P_star_TMS)
% Same as *_Optimized but builds an initial ODE state (no warm-start yet).

    Pprop_initial = 16242000;
    load_percent_ratio = Pprop / Pprop_initial;
    
    eff_FC = [0.484 0.55 0.572 0.646];
    load_percent = [1, 0.7, 0.55, 0.30];
    eff_fc_fit = polyfit(load_percent, eff_FC, 2);
    eff_DCDC = [0.98 0.981 0.982 0.982];
    eff_DCDC_fit = polyfit(load_percent, eff_DCDC, 2);
    eff_HTS = [0.999 0.999 0.999 0.999];
    eff_HTS_fit = polyfit(load_percent, eff_HTS, 2); 
    eff_motorD = [0.99 0.991 0.992 0.992];
    eff_motorD_fit = polyfit(load_percent, eff_motorD, 2);
    eff_motor = [0.992 0.994 0.995 0.997];
    eff_motor_fit = polyfit(load_percent, eff_motor, 2);
    
    v.eff_FC = polyval(eff_fc_fit, load_percent_ratio);
    v.eff_DCDC = polyval(eff_DCDC_fit, load_percent_ratio);
    v.eff_HTS = polyval(eff_HTS_fit, load_percent_ratio);
    v.eff_DCAC = polyval(eff_motorD_fit, load_percent_ratio);
    v.eff_M = polyval(eff_motor_fit, load_percent_ratio);
    
    components = simPar_IZEA_2_Optimized(v, Pprop);
    P_M=components(1);
    Q_dot_M=components(2);
    P_DCAC=components(3);
    Q_dot_DCAC=components(4);
    P_HTS=components(5);
    Q_dot_HTS=components(6);
   
    target_fun = @(x) [
        x(1) - ((P_HTS + x(2) + PEMFCj(x(1), altitude, speed, A_FC_eff)) * ((1 / v.eff_DCDC) - 1)) ...
              - (P_HTS+ x(2) + PEMFCj(x(1), altitude, speed, A_FC_eff));
        x(2) - 0.04 * 0.95 * ((x(1) / v.eff_FC) - x(1))
    ];

    x0 = [P_FC_el_initial, P_star_TMS];

    % fsolve with guard & retry as above
    opts = optimoptions('fsolve', ...
        'Display','off', ...
        'Algorithm','trust-region-dogleg', ...
        'StepTolerance',1e-9, ...
        'FunctionTolerance',1e-9, ...
        'MaxIterations',400, ...
        'MaxFunctionEvaluations',4000);

    [x_solution, fval, exitflag, output] = fsolve(target_fun, x0, opts);

    if exitflag <= 0 || any(~isfinite(x_solution))
        warning(['PowerModel2j_IZEA_2_Initial: fsolve did not converge. ', ...
                 'exitflag=%d, iters=%d, normF=%g. Retrying with Levenberg–Marquardt.'], ...
                 exitflag, output.iterations, norm(fval));

        optsLM = optimoptions(opts, 'Algorithm','levenberg-marquardt');
        [x_solution2, fval2, exitflag2, output2] = fsolve(target_fun, x0, optsLM);

        if exitflag2 > 0 && all(isfinite(x_solution2))
            x_solution = x_solution2;  fval = fval2;  exitflag = exitflag2;  output = output2;
        else
            warning(['PowerModel2j_IZEA_2_Initial: retry also failed. ', ...
                     'Keeping initial guess to continue simulation.']);
            x_solution = x0;  % fallback
        end
    end

    P_FC_el = x_solution(1);
    P_TMS = x_solution(2);

    P_comp= PEMFCj_Optimized(P_FC_el, altitude, speed, A_FC_eff);
    Q_dot_DCDC = (P_HTS + P_comp+ P_TMS) * ((1/v.eff_DCDC) - 1);
    P_DCDC = P_FC_el - Q_dot_DCDC;

    [temperatures,m_dot,Q_dot_FC_H2, yinitial,Re_HTS] = TempCalcj_IZEA_2_Initial(Pprop, P_FC_el, A_FC_eff, Q_dot_DCDC,v);
    P_H2_LHV_eff = P_FC_el ./ v.eff_FC;
    P = [P_M, P_DCAC, P_HTS,P_DCDC, P_FC_el, P_comp,P_TMS, P_H2_LHV_eff];
    Q_dot = [Q_dot_M, Q_dot_DCAC, Q_dot_HTS, Q_dot_DCDC, Q_dot_FC_H2];
    T = temperatures;
    jcell = jcalculation(P_FC_el,A_FC_eff);
end

function components = simPar_IZEA_2_Optimized(v, Pprop)
% Propagate power and compute losses through Motor → DC/AC → HTS

P_Downstr = Pprop;

% Motor stage
Q_dot_M = P_Downstr * ((1./v.eff_M)-1);
P_M = P_Downstr + Q_dot_M;

P_Downstr = P_M;

% DC/AC stage (includes PF on input)
PF_M = 0.89;
Q_dot_DCAC = (P_Downstr./PF_M) * ((1./v.eff_DCAC)-1);
P_DCAC = P_Downstr + Q_dot_DCAC;

P_Downstr = P_DCAC;

% HTS transport
Q_dot_HTS = P_Downstr * ((1./v.eff_HTS)-1);
P_HTS = P_Downstr + Q_dot_HTS;

components=[P_M  Q_dot_M  P_DCAC  Q_dot_DCAC  P_HTS  Q_dot_HTS];
end

function Pcomp = PEMFCj_Optimized(P_FC_el, altitude, speed, AFC)
% Estimate compressor electric power for FC air supply given
% ambient (altitude, speed), target FC power, and area.

    % Constants and ambient model
    ss = 343;  cp_air = 1000;  y_air = 1.4;
    Ma = speed / ss;
    M_air = 0.02896;  F = 96485;  x_O2 = 0.21;  lamda_O2 = 1.8;
    eff_comp_m = 0.97; eff_comp_el = 0.94; eff_comp_pc = 0.95; eff_comp_s = 0.76;
    p_FCHX_in = 175000;  R_sp_air = 287;  eff_pr = 0.75;

    % Current density to get air stoichiometry
    jcell = jcalculation(P_FC_el,AFC);
    
    pamb = 101325 * (1 - 0.0065* (altitude./288.15))^5.2561;
    Tamb = 288.15-6.5*(altitude./1000)+24*(1-(altitude./12664)); % simple lapse + adjustment
    
    % Air mass flow (stoichiometry)
    m_dot_air_in = (jcell * AFC * lamda_O2 * M_air) / (4 * F * x_O2);
    
    % Total inlet, pressure ratio, specific work
    T_comp_in = Tamb + ((sqrt(y_air * R_sp_air * Tamb) * Ma)^2) / (2 * cp_air);
    p_comp_in_stat = pamb * (T_comp_in / Tamb)^(y_air / (y_air - 1));
    p_comp_in = eff_pr * (p_comp_in_stat - pamb) + pamb;
    delta_h_comp = (1 / eff_comp_s) * cp_air * T_comp_in * ((p_FCHX_in / p_comp_in)^(R_sp_air / cp_air) - 1);
    
    % Electrical power with chain efficiencies
    Pcomp = (delta_h_comp * m_dot_air_in) / (eff_comp_m * eff_comp_el * eff_comp_pc);
end

function [temperatures,m_dot,Q_dot_FC_H2, y,Re_HTS] = TempCalcj_IZEA_2_Optimized(Pprop, P_FC_el, AFC, Q_dot_DCDC, v,yinitial)
% Transient thermal model for 4 serial HXs with H2/He channels and wall.
% Warm-started from previous ODE state 'yinitial'.

    % Fluids/pressures/geometry
    p.ele1 = 'H2';
    p.ele2 = 'He';
    p.p_ele1 = 166000;                 % [Pa] H2
    p.p_ele2 = 4e05*ones(1,5);         % [Pa] He (per component)
    p.rhow = 8960;                     % [kg/m^3] wall density (Cu)
    % p.cpw, p.kw computed via functions
    p.thx = 0.002;   p.hc_ele1 = 0.005;  p.hc_ele2 = 0.01;
    p.H = 10;        p.W = 1;            p.nx = 5;

    % Power-chain heat loads by block
    components = simPar_IZEA_2_Optimized(v, Pprop);
    Q_dot_M = components(2);
    Q_dot_DCAC = components(4);
    Q_dot_HTS = components(6);
  
    % Cathode/stoichiometry
    F = 96485; x_O2 = 0.21; lamda_H2_net = 1.05; lamda_O2 = 1.8; M_air = 0.02896;
    jcell = jcalculation(P_FC_el,AFC);

    % Mass flows
    Fc_eff = v.eff_FC;  LHV_H2 = 120E6;
    p.mdot_ele1 = (lamda_H2_net * P_FC_el) / (Fc_eff * LHV_H2);        % H2
    m_dot_air_in = (jcell * AFC * lamda_O2 * M_air) / (4 * F * x_O2);  % Air
     
    % MOL setup for this mission step
    tf = 30;  x = linspace(0, p.H, p.nx);  dx = diff(x);  p.dx = dx(1);  tspan = [0 tf];
    options = odeset('Reltol', 1e-2, 'AbsTol', 1e-3);
    
    % Heat loads per component order: [HTS, DC/DC, DC/AC, M]
    p.Q_dot = [Q_dot_HTS Q_dot_DCDC Q_dot_DCAC Q_dot_M];
    
    % Warm-start temperatures from last state
    p.T_H2_in = [yinitial(end,1) yinitial(end,1+p.nx) yinitial(end,1+2*p.nx) yinitial(end,1+3*p.nx)];
    T_H2 = yinitial(end,1:4*p.nx);
    T_He = yinitial(end,4*p.nx+1:8*p.nx);
    T_Wall = yinitial(end,8*p.nx+1:12*p.nx);
    yc = [T_H2'; T_He'; T_Wall'];

    % Integrate
    [~, y] = ode15s(@(t, y) ode_optimized(t, y, p), tspan, yc, options);

    % Mass flows out
    m_dot_H2_in = p.mdot_ele1;
    m_dot = [m_dot_H2_in, m_dot_air_in];

    % Outlet temps (last axial node of each component)
    T_HTS = y(:, p.nx);
    T_DCDC = y(:, 2*p.nx);
    T_DCAC = y(:, 3*p.nx);
    T_M = y(:, 4*p.nx);

    temperatures = [T_HTS, T_DCDC, T_DCAC, T_M];

    % Reynolds at HTS outlet (consistent with diagnostic label)
    A=p.hc_ele1 * p.W;
    dh=(2*p.hc_ele1 * p.W)./(p.hc_ele1 + p.W);
    mu_HTS=py.CoolProp.CoolProp.PropsSI('V', 'P', 166000, 'T', T_HTS(end), p.ele1);
    Re_HTS = (m_dot_H2_in * dh)./(A* mu_HTS);

    % Residual FC heat absorbed by H2 en route to PEM
    C_h2 = CalcCap(T_M(end));  T_PEM = 358; % [K]
    Q_dot_FC_H2 = C_h2 * m_dot_H2_in * (T_PEM - T_M(end));
end

function [temperatures,m_dot, Q_dot_FC_H2,yinitial,Re_HTS] = TempCalcj_IZEA_2_Initial(Pprop, P_FC_el, AFC, Q_dot_DCDC, v)
% Same as above but builds an initial uniform state (no warm-start yet).

    % Fluids/pressures/geometry
    p.ele1 = 'H2';
    p.ele2 = 'He';
    p.p_ele1 = 166000;
    p.p_ele2 = 4e5*ones(1,5);
    p.rhow = 8960;
    % p.cpw, p.kw via functions
    p.h_ele1_w = 150;  % (not used later; retained)
    p.thx = 0.002;  p.hc_ele1 = 0.005;  p.hc_ele2 = 0.01;
    p.H = 10;  p.W = 1;  p.nx = 5;
    
    components = simPar_IZEA_2_Optimized(v, Pprop);
    Q_dot_M = components(2);
    Q_dot_DCAC = components(4);
    Q_dot_HTS = components(6);
    
    % Stoichiometry
    F = 96485;  x_O2 = 0.21;  lamda_H2_net = 1.05;  lamda_O2 = 1.8;  M_air = 0.02896;
    jcell = jcalculation(P_FC_el,AFC);
    
    % Mass flows
    Fc_eff = v.eff_FC;  LHV_H2 = 120E6;
    p.mdot_ele1 = (lamda_H2_net * P_FC_el) / (Fc_eff * LHV_H2);
    m_dot_air_in = (jcell * AFC * lamda_O2 * M_air) / (4 * F * x_O2);
    
    % MOL setup
    tf = 30;  x = linspace(0, p.H, p.nx);  dx = diff(x); p.dx = dx(1); tspan = [0 tf];
    options = odeset('Reltol', 1e-2, 'AbsTol', 1e-3);
    
    % Heat loads [HTS, DC/DC, DC/AC, M]
    p.Q_dot = [Q_dot_HTS Q_dot_DCDC Q_dot_DCAC Q_dot_M];

    % Initial uniform temps (H2=23 K, He=30 K, Wall average)
    T_H2=23*ones(4*p.nx,1);
    T_He=[30;30;30;30;30;30;30;30;30;30;30;30;30;30;30;30;30;30;30;30];
    T_Wall=(T_He+T_H2)./2;
    y_initial = [T_H2; T_He; T_Wall];

    [~, yinitial] = ode15s(@(t, y) ode_optimized(t, y, p), tspan, y_initial, options);
    
    % Mass flows out
    m_dot_H2_in = p.mdot_ele1;
    m_dot = [m_dot_H2_in, m_dot_air_in];

    % Outlet temps
    T_HTS = yinitial(:, p.nx);
    T_DCDC = yinitial(:, 2*p.nx);
    T_DCAC = yinitial(:, 3*p.nx);
    T_M = yinitial(:, 4*p.nx);
    temperatures = [T_HTS, T_DCDC, T_DCAC, T_M];

    % Reynolds at HTS outlet
    A=p.hc_ele1 * p.W;
    dh=(2*p.hc_ele1 * p.W)./(p.hc_ele1 + p.W);
    mu_HTS=py.CoolProp.CoolProp.PropsSI('V', 'P', 166000, 'T', T_HTS(end), p.ele1);
    Re_HTS = (m_dot_H2_in * dh)./(A* mu_HTS);

    % Residual FC heat in H2
    C_h2 = CalcCap(T_M(end));  T_PEM = 358;
    Q_dot_FC_H2 = C_h2 * m_dot_H2_in * (T_PEM - T_M(end));
end

function Temp = ode_optimized(t, y, p)
% Semi-discrete energy balances for 4 components × (H2, He, Wall).
% - Convection in H2/He channels + wall conduction + convective coupling
% - Pressure drop marching to update density/velocity
% - Gnielinski/Churchill correlations, with laminar/transition/turbulent logic

    n_comp = 4;

    % Derivative buffers
    dT_ele1dx = zeros(p.nx, 1);
    dT_ele2dx = zeros(p.nx, 1);
    dT2_walldx2 = zeros(p.nx, 1); 
    dT_ele1dt = zeros(n_comp*p.nx,1);
    dT_ele2dt = zeros(n_comp*p.nx,1);
    dTwdt = zeros(n_comp*p.nx,1);

    % Unpack state
    T_H2 = y(1:p.nx*n_comp);
    T_He = y(n_comp*p.nx+1:2*n_comp*p.nx);
    T_Wall = y(2*n_comp*p.nx+1:3*n_comp*p.nx);

    % Simple CoolProp cache (reduces Python call overhead)
    persistent coolprop_cache
    if isempty(coolprop_cache)
        coolprop_cache = containers.Map('KeyType', 'char', 'ValueType', 'double');
    end
    function val = getCoolPropValue(prop, P, T, fluid)
        key = sprintf('%s_%f_%f_%s', prop, P, T, fluid);
        if ~isKey(coolprop_cache, key)
            val = double(py.CoolProp.CoolProp.PropsSI(prop, 'P', P, 'T', T, fluid));
            if isnan(val) || isinf(val)
                error("CoolProp error with P = %.2f, T = %.2f", P, T);
            end
            coolprop_cache(key) = val;
        else
            val = coolprop_cache(key);
        end
    end

    % Nested helper: H2 inlet into component k (serial chaining)
    function T_in_k = get_H2_inlet(k)
        if k == 1
            T_in_k = 23;           % first HX inlet (fixed guess)
        else
            T_in_k = T_H2((k-1)*p.nx); % last node of previous HX
        end
    end

    T_max=[60,140,140,140];  % He inlet cap per HX [K] (control/limit)
    %T_min=[50,110,110,110]; % (unused)
    dP_H2_total=zeros(n_comp,1);
    dP_He_total=zeros(n_comp,1);

    for k = 1:n_comp
        % Axial indices for this component
        i_start = 1 + (k-1)*p.nx;
        i_end = k*p.nx;

        % H2 inlet temperature for this HX
        T_in_k = get_H2_inlet(k);

        % Allocate property arrays (H2 + wall)
        p.rho_ele1=zeros(p.nx,1);
        p.cp_ele1=zeros(p.nx,1);
        p.mu_ele1=zeros(p.nx,1);
        p.k_ele1=zeros(p.nx,1);
        p.u_ele1=zeros(p.nx,1);
        p.cpw =zeros(p.nx,1);  % wall heat cap
        p.kw =zeros(p.nx,1);   % wall conductivity

        % Fill H2/wall properties
        for i = 1:p.nx
            p.rho_ele1(i) = getCoolPropValue('D', p.p_ele1,T_H2(i+(k-1)*p.nx), p.ele1);
            p.cp_ele1(i) = getCoolPropValue('C', p.p_ele1, T_H2(i+(k-1)*p.nx), p.ele1);
            p.mu_ele1(i) = getCoolPropValue('V', p.p_ele1, T_H2(i+(k-1)*p.nx), p.ele1);
            p.k_ele1(i) = getCoolPropValue('L', p.p_ele1, T_H2(i+(k-1)*p.nx), p.ele1);
            p.cpw(i)=Wall_capacity(T_Wall(i+(k-1)*p.nx));
            p.kw(i)=Wall_conductivity(T_Wall(i+(k-1)*p.nx));
        end

        % H2 velocity from mass conservation
        p.u_ele1 = p.mdot_ele1 ./ (p.rho_ele1 .* p.hc_ele1.*p.W);

        % Average properties for correlations
        rho_H2_avg=mean(p.rho_ele1);
        cp_H2_avg=mean(p.cp_ele1);
        mu_H2_avg=mean(p.mu_ele1);
        k_H2_avg=mean(p.k_ele1);
        u_H2_avg=mean(p.u_ele1);
        
        % Hydraulic diameter and Reynolds (H2)
        epsilon= 1.5e-05;
        Dh = (2*p.hc_ele1*p.W)/(p.hc_ele1+p.W);
        ReD = (u_H2_avg * Dh * rho_H2_avg) / mu_H2_avg;

        % Friction factor f (H2)
        if ReD < 2300
            f= 96 / ReD;
        elseif ReD >= 2300 && ReD <= 10000
            A = (2.457*log( (7/ReD)^0.9 + 0.27*(epsilon./Dh) ) )^16;
            B = (37530/ReD)^16;
            f = 8*((8/ReD)^12 + 1/(A+B)^(1.5))^(1/12); 
        else 
            f=(1.8*log10(ReD)-1.5).^-2;
        end

        % Pressure drop marching (H2)
        dP_H2=zeros(p.nx,1);
        localp_H2=p.p_ele1;
        for i = 1:(p.nx-1)
            dP_H2(i)= 0.5*p.rho_ele1(i)*(p.u_ele1(i).^2)*(f*p.dx./Dh + 2*(p.u_ele1(i+1)-p.u_ele1(i))/p.u_ele1(i));
            localp_H2 = localp_H2 - dP_H2(i);
            p.rho_ele1(i+1) = getCoolPropValue('D', localp_H2,T_H2(i+(k-1)*p.nx), p.ele1);
            p.u_ele1(i+1) = p.mdot_ele1 ./ (p.rho_ele1(i+1) .* p.hc_ele1.*p.W);
        end
        dP_H2_total(k)=sum(dP_H2);
        
        % Nusselt (H2) by regime
        Pr = cp_H2_avg * mu_H2_avg / k_H2_avg;
        if ReD <= 2300
            fprintf("Laminar Region H2");
            Nu=7.541; % laminar surrogate for short L/D
        elseif ReD>= 2300 && ReD<=10000
            fprintf("Transition Region H2");
            l=(ReD-2300)./(10000-2300);
            Nu_lam_2300=7.541;
            Nu_turb_10000=((0.0308/8)*10000*Pr)/(1+12.7*sqrt(0.0308/8)*Pr.^(2/3)-1);
            Nu=(1-l)*Nu_lam_2300 + l*Nu_turb_10000;
        else
            Nu = ((f/8) * (ReD) * Pr) / (1 + 12.7*(f/8)^0.5*(Pr^(2/3) - 1))*(mean(T_H2)./mean(T_Wall)).^-0.2;
        end
        h_ele1 = (Nu * k_H2_avg) / Dh;

        % Axial derivative (upwind at inlet)
        dT_ele1dx(2:end) = (T_H2(i_start+1:i_end) - T_H2(i_start:i_end-1)) / p.dx;
        if k == 1
            dT_ele1dx(1) = (T_H2(i_start) - T_in_k) / (p.dx/2);
        else
            dT_ele1dx(1) = (T_H2(i_start) - T_H2(i_start-1)) / p.dx;
        end

        % H2 energy equation
        dT_ele1dt(i_start:i_end) = - p.u_ele1 .* dT_ele1dx + (h_ele1/p.hc_ele1) * (T_Wall(i_start:i_end) - T_H2(i_start:i_end))./ (p.rho_ele1 .* p.cp_ele1);

        % He side inlet cap/target
        T_He_in_k = T_max(k);

        % He properties arrays
        p.rho_ele2=zeros(p.nx,1);
        p.cp_ele2=zeros(p.nx,1);
        p.mu_ele2=zeros(p.nx,1);
        p.k_ele2=zeros(p.nx,1);
        p.u_ele2=zeros(p.nx,1);

        % He mass flow to meet heat load (with cap at 1 kg/s)
        h_out_ele2 = py.CoolProp.CoolProp.PropsSI('H', 'P',p.p_ele2(k),'T', T_He(i_start), p.ele2);
        h_in_ele2 = py.CoolProp.CoolProp.PropsSI('H', 'P',p.p_ele2(k),'T', T_max(k), p.ele2);
        delta_h = h_in_ele2 - h_out_ele2;
        p.mdot_ele2(k) = p.Q_dot(k) / abs(delta_h);
        if p.mdot_ele2(k) > 1
            p.mdot_ele2(k)=1;
            h_in_ele2 = h_out_ele2 + (p.Q_dot(k)/p.mdot_ele2(k));
            T_He_in_k= py.CoolProp.CoolProp.PropsSI('T', 'P', p.p_ele2(k), 'H', h_in_ele2, p.ele2);
        end

        % Fill He properties
        for i = 1:p.nx
            p.rho_ele2(i) = getCoolPropValue('D', p.p_ele2(k), T_He(i+(k-1)*p.nx), p.ele2);
            p.cp_ele2(i) = getCoolPropValue('C', p.p_ele2(k), T_He(i+(k-1)*p.nx), p.ele2);
            p.mu_ele2(i) = getCoolPropValue('V', p.p_ele2(k), T_He(i+(k-1)*p.nx), p.ele2);
            p.k_ele2(i) = getCoolPropValue('L', p.p_ele2(k), T_He(i+(k-1)*p.nx), p.ele2);
        end
        p.u_ele2 = p.mdot_ele2(k) ./ (p.rho_ele2 .* p.hc_ele2.*p.W);

        % Average He properties
        rho_He_avg=mean(p.rho_ele2);
        cp_He_avg=mean(p.cp_ele2);
        mu_He_avg=mean(p.mu_ele2);
        k_He_avg=mean(p.k_ele2);
        u_He_avg=mean(p.u_ele2);
        
        % He hydraulic diameter and Re
        Dh_He = (2*p.hc_ele2*p.W)/(p.hc_ele2+p.W);
        ReD2 = (u_He_avg * Dh_He * rho_He_avg) / mu_He_avg;

        % Friction factor f2 (He)
        if ReD2 < 2300
            f2 = 96 / ReD2;
        elseif ReD2 >=2300 && ReD2 <=10000
            A = (2.457*log( (7/ReD2)^0.9 + 0.27*(epsilon./Dh_He) ) )^16;
            B = (37530/ReD2)^16;
            f2 = 8*((8/ReD2)^12 + 1/(A+B)^(1.5))^(1/12);
        else 
            f2=(1.8*log10(ReD2)-1.5).^-2;
        end

        % He pressure drop marching
        dP_He=zeros(p.nx,1);
        localp_He=zeros(n_comp,1);
        localp_He(k) = p.p_ele2(k);
        for i = 1:(p.nx-1)
            dP_He(i)= 0.5*p.rho_ele2(i)*p.u_ele2(i).^2*(f2*p.dx./Dh_He + 2*(p.u_ele2(i+1)-p.u_ele2(i))/p.u_ele2(i));
            localp_He(k) = localp_He(k) - dP_He(i);
            p.rho_ele2(i+1) = getCoolPropValue('D', localp_He(k), T_He(i+(k-1)*p.nx), p.ele2);
            p.u_ele2(i+1) = p.mdot_ele2(k) ./ (p.rho_ele2(i+1) .* p.hc_ele2.*p.W);
        end
        dP_He_total(k)=sum(dP_He);

        if any(~isfinite([p.rho_ele1;p.rho_ele2]))
            error('CoolProp returned NAN—check P/T range')
        end
       
        % He Nusselt by regime
        Pr2 = cp_He_avg * mu_He_avg / k_He_avg;
        if ReD2 <= 2300
            fprintf("Laminar Region He");
            NuD2=7.541;
        elseif ReD2>= 2300 && ReD2<=10000
            fprintf("Transition Region He");
            l=(ReD2-2300)./(10000-2300);
            Nu_lam_2300=7.541;
            Nu_turb_10000=((0.0308/8)*10000*Pr2)/(1+12.7*sqrt(0.0308/8)*Pr2.^(2/3)-1);
            NuD2=(1-l)*Nu_lam_2300 + l*Nu_turb_10000;
        else
            NuD2 = ((f2/8) * (ReD2) * Pr2) / (1 + 12.7*(f2/8)^0.5*(Pr2^(2/3) - 1));
        end
        h_ele2 = (NuD2 * k_He_avg) / Dh_He;

        % Backward differencing (outlet boundary at i_end)
        dT_ele2dx(1:end-1) = (T_He(i_start:i_end-1) - T_He(i_start+1:i_end)) / p.dx;
        dT_ele2dx(end) = (T_He(i_end) - T_He_in_k) / (p.dx/2);

        % He energy equation
        dT_ele2dt(i_start:i_end) =  -p.u_ele2 .* dT_ele2dx + (h_ele2/p.hc_ele2) * (T_Wall(i_start:i_end) - T_He(i_start:i_end))./ (p.rho_ele2 .* p.cp_ele2);
            
        % Wall conduction (second derivative) with one-sided ends
        dT2_walldx2(2:p.nx-1) = (T_Wall(i_start+2:i_end) + T_Wall(i_start:i_end-2) ...
                - 2*T_Wall(i_start+1:i_end-1)) / (p.dx^2);
        dT2_walldx2(1) = (   2*T_Wall(i_start) - 5*T_Wall(i_start+1) + 4*T_Wall(i_start+2) - T_Wall(i_start+3) ) /(p.dx)^2;
        dT2_walldx2(end) = ( 2*T_Wall(i_end) - 5*T_Wall(i_end-1) + 4*T_Wall(i_end-2) - T_Wall(i_end-3) ) /(p.dx)^2;

        % Wall energy equation (conduction + two-side convection)
        dTwdt(i_start:i_end) = (p.kw .* dT2_walldx2 ...
            + (h_ele1/p.hc_ele1) * (T_H2(i_start:i_end) - T_Wall(i_start:i_end))  + ...
            (h_ele2/p.hc_ele2) * (T_He(i_start:i_end) - T_Wall(i_start:i_end)) )./ (p.rhow .* p.cpw);
    end

    % Return concatenated time derivatives
    Temp = [dT_ele1dt; dT_ele2dt; dTwdt];
end

function C = CalcCap(T)
% Hydrogen cp(T) via CoolProp at 166 kPa
P= 166000 ; % [Pa]
C = py.CoolProp.CoolProp.PropsSI('C','T',T,'P',P,'Hydrogen');
end

function Cw = Wall_capacity(T)
% Empirical wall heat capacity [J/kg-K] vs T
Cw = 10.^(1.131-9.454*(log10(T))+12.99*(log10(T)).^2 -5.501*(log10(T)).^3 + 0.7637*(log10(T)).^4);
end

function kw = Wall_conductivity(T)
% Empirical wall thermal conductivity [W/m-K] vs T
a = 1.97365e4;    % W/mK
b = 6.55194e-2;   % 1/K
c = 4.18166e2;    % W/mK
kw=a.*exp(-b.*T) + c;
end
