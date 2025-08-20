%% ========================================================================
%  Thesis Sim Code — IZEA (PEMFC + HTS + TMS)
%  Main sweep over maximum cell current density (j_cell,max) to find the
%  architecture with minimum stacked system mass. Generates plots and
%  exports a figure.
%
%  INPUT DATA (Excel):
%    Mission_Profile_Data_Total.xlsx with columns:
%      - time_sec_       [s]
%      - altitude_m_     [m]
%      - speed_ms_       [m/s]
%      - P_W_            [W]  <-- propulsive power demand
%
%  OUTPUTS:
%    - Figure 1: Stacked component masses vs j_cell,max
%    - min_mass_onlyFC.jpg (300 dpi)
%    - Figure 2: Mission current-density profile
%
%  NOTE:
%    This header and all added lines are comments/documentation only.
%    No executable code has been altered.
%  ========================================================================

%% Import DATA values
clear; clc;
data = readtable('Mission_Profile_Data_Total.xlsx'); % read the mission profile
time = data.time_sec_;
% time_seg = data.timeSegment_sec_;
time_step = length(time);
altitude = data.altitude_m_;
speed= data.speed_ms_;
P_prop = data.P_W_;
P_prop_initial = P_prop(1);

% Extract first time-step operating point (for illustrative initialization)
t=1;
Pprop = P_prop(t);   % [W] propulsive power at t=1
h = altitude(t);     % [m]
s = speed(t);        % [m/s]    % obtain altitude, speed, and power for this time

% Sweep bounds for maximum cell current density j_cell,max [A/cm^2]
jmin = 0.1;
jmax = 1.6;
nstep = 0.1; % step for j_cell,max

%% Set empty matrices

range = jmin:nstep:jmax;            % vector of j_cell,max sweep points
n_iter = length(range);             % number of sweep points

% Stacked masses per component across the sweep (cumulative layering)
all_masses = zeros(9,n_iter);
line_zero = zeros(1,n_iter);

% Per-component absolute mass records (for convenience/plotting)
m_stack_all = zeros(1,n_iter);
m_comp_all  = zeros(1,n_iter);
m_hum_all   = zeros(1,n_iter);
m_BoP_all   = zeros(1,n_iter);
m_TMS_all   = zeros(1,n_iter);
m_DCDC_all  = zeros(1,n_iter);
m_HTS_all   = zeros(1,n_iter);
m_DCAC_all  = zeros(1,n_iter);
m_SCM_all   = zeros(1,n_iter);

% Fuel-cell effective area for each sweep point
all_A_FC = zeros(n_iter,1);   % stores A_FC_eff [cm^2] (effective active area)

%% Iteration process to find optimal current density
% Fixed reference electrical power levels used inside the models (W)
P_FC=1.8e07;
P_star_TMS=5e06;

for i = 1:n_iter

    % --------------------------------------------------------------------
    % 1) Given j_cell,max, altitude, speed, and P_prop, compute power splits,
    %    heat loads, air mass flow, and required FC active area A_FC_eff.
    % --------------------------------------------------------------------
    [P,Q_dot,m_dot,A_FC_eff] = PowerModel2_IZEA_light(Pprop,h,s,jmin + nstep*(i-1),P_FC,P_star_TMS); 
    % NOTE: PowerModel2_* version takes j_cell,max as input.

    % --------------------------------------------------------------------
    % 2) Run full-mission iteration with that A_FC_eff to get time histories
    %    of powers/heat and the maxima needed for sizing and mass calc.
    % --------------------------------------------------------------------
    [Powers,Qtransfer,j,Maximum] = MissionIteration_IZEA_light(A_FC_eff,P_prop);
    % Maximum contains peak ratings required for component mass sizing.

    %% These maximum powers are needed in the massCalc function
    P_M            = Maximum(1); % Motor electrical power rating [W]
    P_DCAC         = Maximum(2); % Inverter (DC/AC) rating [W]
    P_HTS          = Maximum(3); % HTS cable through-power [W]
    P_DCDC         = Maximum(4); % DC/DC rating [W]
    P_FC_LHV_eff   = Maximum(5); % Effective H2 LHV power into FC [W] (used as proxy)
    P_comp         = Maximum(6); % Compressor power [W]
    P_TMS          = Maximum(7); % Thermal management system power [W]
    m_dot_air      = Maximum(8); % Air mass flow into FC cathode [kg/s]

    %% Calculate overall mass of the components (cumulative stack)
    [m_stack,m_comp,m_hum,m_BoP,m_TMS,m_DCDC,m_HTS,m_DCAC,m_SCM] = ...
        massCalc(P_M,P_HTS,P_DCAC,P_DCDC,A_FC_eff,P_TMS,P_comp,m_dot_air); % component masses [kg]

    % Store absolute masses
    m_stack_all (i) =  m_stack;
    m_comp_all  (i) =  m_comp;
    m_hum_all   (i) =  m_hum;
    m_BoP_all   (i) =  m_BoP;
    m_TMS_all   (i) =  m_TMS;
    m_DCDC_all  (i) =  m_DCDC;
    m_HTS_all   (i) =  m_HTS;
    m_DCAC_all  (i) =  m_DCAC;
    m_SCM_all   (i) =  m_SCM;

    % Store cumulative layering for stacked-area plotting
    all_masses(:,i) = [m_stack;m_comp;m_hum;m_BoP;m_TMS;m_DCDC;m_HTS;m_DCAC;m_SCM];

    % Store effective FC area
    all_A_FC(i) = A_FC_eff;
end
 
%% Calculates minimum mass across the sweep

all_masses = [line_zero;all_masses];        % prepend zero baseline for fill
all_total_masses = all_masses(9,:);         % total mass trace for plotting
[prop_sys_min_mass,idx] = min(all_total_masses);
jcell_max = 0.2 + 0.1*idx;                  % derived j_cell,max from index (legacy)

% NOTE:
% all_masses is (N_components+1 x N_jcell) = 10 x n_iter here
% range is the array of j_cell,max values

%% ===== Figure 1: Stacked component masses vs j_cell,max =================
figure(1); clf
hold on

% Component labels (stack order must match all_masses layering)
component_labels = {
    'PEMFC stack', ...
    'PEMFC Humidifier', ...
    'PEMFC compressor', ...
    'PEMFC BoP', ...
    'TMS', ...
    'DC/DC converter', ...
    'HTS cable', ...
    'DC/AC inverter', ...
    'Motor'
};

% Color palette for stacked fill (RGB 0-1)
component_colors = [
    0.85, 0.33, 0.30;   % PEMFC stack
    0.94, 0.67, 0.66;   % Humidifier
    0.80, 0.39, 0.45;   % Compressor
    0.6,  0.2,  0.3;    % BoP
    0.5,  0.7,  1.0;    % TMS
    0.1,  0.4,  0.8;    % DC/DC
    1.0,  0.4,  0.0;    % HTS
    0.77, 0.8,  0.9;    % DC/AC
    0.2,  0.6,  0.2;    % Motor
];

% Draw stacked filled layers between component curves
for k = 1:size(all_masses, 1)-1
    X = [range, fliplr(range)];
    Y1 = all_masses(k, :);
    Y2 = all_masses(k+1, :);
    Y = [Y1, fliplr(Y2)];
    fill(X, Y, component_colors(k,:), 'EdgeColor', 'none', 'FaceAlpha', 0.85);
end

% Outline of total mass
plot(range, all_masses(end,:), 'k--', 'LineWidth', 2);

% Mark the minimum total mass and vertical guide
[mtot_min, idx_min] = min(all_masses(end,:));
jcell_opt = range(idx_min);
plot(jcell_opt, mtot_min, 'k.', 'MarkerSize', 20);
xline(jcell_opt, 'k--', 'LineWidth', 1.5);
text(jcell_opt + 0.01, mtot_min + 50, ...
    sprintf('m_{tot,min}', round(mtot_min)), ...
    'FontSize', 14, 'FontWeight', 'bold');

% Axes, labels, and formatting
xlabel("Maximum cell current density (j_{cell,max}) [A/cm^2]", ...
    "FontSize", 18, "FontWeight", "bold");
ylabel("Stacked weights [kg]", ...
    "FontSize", 18, "FontWeight", "bold");
set(gca, "FontSize", 18, "FontWeight", "bold");
grid on
set(gca, 'GridAlpha', 0.4, 'GridLineStyle', '-', ...
    'GridColor', [0.3 0.3 0.3], 'GridLineWidth', 0.75);

% Legend outside for readability
legend(component_labels, 'Location', 'eastoutside', ...
    'FontSize', 18, 'Box', 'off');

% Tight layout and export
xlim([min(range), max(range)]);
ylim([0, max(all_masses(end,:)) + 300]);
title('Component Masses vs. Maximum Cell Current Density', 'FontSize', 18);
% print('min_mass_onlyFC.jpg','-djpeg');%,'-r300'
exportgraphics(gcf,'min_mass_onlyFC.jpg','Resolution',300);

%% ===== Figure 2: Mission current-density profile ========================
figure(2)
plot(time/3600,j,LineWidth=2,Color=[0 1 0]);   
hold on
set(gca, 'LineWidth', 2);
title("Current density","FontSize",20,"FontWeight","bold");
% xlabel("Time [s]");
xticks([2 4 6 8]);
ylabel("j[A/cm2]","FontSize",15,"FontWeight","bold");
set(gca,"FontSize",15,"FontWeight","bold")
grid on
grid(gca, 'on');
set(gca, 'GridAlpha', 0.5);      % grid transparency
set(gca, 'GridLineStyle', '-');  % grid style
set(gca, 'GridColor', 'k');      % grid color
set(gca, 'GridLineWidth', 1);    % grid width
hold off


%% ========================================================================
%  MissionIteration_IZEA_light
%  ------------------------------------------------------------------------
%  PURPOSE:
%    Runs the mission profile using a fixed FC effective area (A_FC_eff),
%    computing power flows, heat transfers, current density, and tracking
%    the maximum ratings needed for component sizing.
%
%  INPUTS:
%    A_FC_eff  [cm^2]   - Effective PEMFC active area
%    P_prop    [W]      - Vector of mission propulsive power demand
%
%  OUTPUTS:
%    Powers    [n x 8]  - Per-timestep power channels
%    Qtransfer [n x 5]  - Per-timestep heat transfer channels
%    j         [n x 1]  - Cell current density time series [A/cm^2]
%    Maximum   [1 x 8]  - Peak ratings across mission:
%                         [P_M, P_DCAC, P_HTS, P_DCDC, P_FC_LHV_eff, ...
%                          P_comp, P_TMS, m_dot_air]
%  ========================================================================
function [Powers,Qtransfer,j,Maximum] = MissionIteration_IZEA_light(A_FC_eff,P_prop)
%% 1 Part we get the data from the excell

data = readtable('Mission_Profile_Data_Total.xlsx');

%% 2 Part we take the parameters and plot them against time

time = data.time_sec_;
altitude = data.altitude_m_;
speed= data.speed_ms_;

n=size(time,1);
Powers = zeros(n,8);
Qtransfer = zeros(n,5);
j = zeros(n,1);

% Initialize maxima trackers for sizing
P_SCM_max= 0;
P_DCAC_max=0;
P_HTS_max=0;
P_DCDC_max= 0 ;
P_FC_LHV_eff_max= 0;
P_comp_max = 0;
P_TMS_max=0;
m_dot_air_max = 0;

%% 3 Part Mission Iteration model

P_FC_el = 1.8e07;     % initial FC electrical power guess [W]
P_star_TMS = 5e06;    % initial TMS power guess [W]

for i=1:(n-1)
    % t = time(i);
    h = altitude(i);
    s = speed(i);
    p = P_prop(i);

    % Compute point solution with j as output (A_FC_eff fixed)
    [P,Q_dot,m_dot,jcell] = PowerModel2j_IZEA_light(p,h,s,A_FC_eff,P_FC_el,P_star_TMS);

    % Unpack channels for clarity
    P_SCM = P(1);
    P_DCAC = P(2);
    P_HTS = P(3);
    P_DCDC = P(4);
    P_FC_el = P(5);
    P_comp = P(6);
    P_TMS  = P(7);
    P_FC_LHV_eff = P(8);
    m_dot_air = m_dot;

    % Use current TMS power as next-iteration seed
    P_star_TMS=P_TMS;

    % Track maxima
    if P_SCM > P_SCM_max,          P_SCM_max = P_SCM;             end
    if P_DCAC > P_DCAC_max,        P_DCAC_max = P_DCAC;           end
    if P_HTS  > P_HTS_max,         P_HTS_max  = P_HTS;            end
    if P_DCDC > P_DCDC_max,        P_DCDC_max = P_DCDC;           end
    if P_FC_LHV_eff > P_FC_LHV_eff_max, P_FC_LHV_eff_max =  P_FC_LHV_eff; end
    if P_comp > P_comp_max,        P_comp_max = P_comp;           end
    if P_TMS  > P_TMS_max,         P_TMS_max  = P_TMS;            end
    if m_dot_air > m_dot_air_max,  m_dot_air_max = m_dot_air;     end

    Powers(i,:)   = P;
    Qtransfer(i,:) = Q_dot;
    j(i) = jcell;

    % Maxima vector
    Maximum = [P_SCM_max P_DCAC_max P_HTS_max P_DCDC_max ...
               P_FC_LHV_eff_max P_comp_max P_TMS_max m_dot_air_max];
end
end

%% ========================================================================
%  PowerModel2_IZEA_light
%  ------------------------------------------------------------------------
%  PURPOSE:
%    Point model that, for a given (Pprop, h, s, jcell_max), solves for the
%    FC electrical power and TMS power consistent with downstream loads and
%    efficiencies. Returns power channels, heat loads, air mass flow, and
%    required FC effective area A_FC_eff.
%
%  INPUTS:
%    Pprop              [W]      Propulsive power demand at this step
%    altitude           [m]
%    speed              [m/s]
%    jcell_max          [A/cm^2] Upper bound for cell current density
%    P_FC_el_initial    [W]      Initial guess for fsolve
%    P_star_TMS         [W]      Initial guess for TMS power
%
%  OUTPUTS:
%    P        [1x7]   = [P_M, P_DCAC, P_HTS, P_DCDC, P_FC_el, P_comp, P_TMS]
%    Q_dot    [1x5]   = [Q_M, Q_DCAC, Q_HTS, Q_DCDC, Q_FC_stack]
%    m_dot    [kg/s]  = Air mass flow into FC cathode
%    A_FC_eff [cm^2]  = Effective FC active area required at jcell_max
%  ========================================================================
function [P,Q_dot,m_dot_air,A_FC_eff] = PowerModel2_IZEA_light(Pprop,altitude,speed,jcell_max,P_FC_el_initial,P_star_TMS)

% Initial Parameters & efficiency maps (quadratic fits vs load fraction)
Pprop_initial = 16242000;

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

    % Evaluate efficiencies at current load fraction
    v.eff_FC   = polyval(eff_fc_fit,   Pprop / Pprop_initial);
    v.eff_DCDC = polyval(eff_DCDC_fit, Pprop / Pprop_initial);
    v.eff_HTS  = polyval(eff_HTS_fit,  Pprop / Pprop_initial);
    v.eff_DCAC = polyval(eff_motorD_fit, Pprop / Pprop_initial);
    v.eff_M    = polyval(eff_motor_fit,  Pprop / Pprop_initial);

% First Calculation: propagate downstream power/heat up to HTS
components   = simPar_IZEA_light(v,Pprop);
P_M          = components(1);
Q_dot_M      = components(2);
P_DCAC       = components(3);
Q_dot_DCAC   = components(4);
P_HTS        = components(5);
Q_dot_HTS    = components(6);

% Solve for FC electrical power and TMS power that balance the chain
x0 = [P_FC_el_initial, P_star_TMS];
x_solution = fsolve(@target_fun, x0);
P_FC_el = x_solution(1);
P_TMS    = x_solution(2);

% Compressor + air mass flow based on jcell_max and flight point
[P_comp, m_dot_air] = PEMFC(P_FC_el, altitude, speed, jcell_max);

% DC/DC losses from HTS + (comp + TMS) branch
Q_dot_DCDC = (P_HTS + (P_comp + P_TMS)) .*((1./v.eff_DCDC)-1);
P_DCDC = P_FC_el - Q_dot_DCDC;

% FC heat and TMS closure (empirical factors)
P_H2_LHV_eff = P_FC_el / v.eff_FC;
Q_dot_FC_stack = 0.95 * (P_H2_LHV_eff - P_FC_el);
P_TMS = 0.04 * Q_dot_FC_stack;

    % Inner function: balances FC output with downstream loads
    function F = target_fun(x)
        [P_comp, ~] = PEMFC(x(1), altitude, speed, jcell_max);

        eq1 = x(1) - ((P_HTS + x(2) + P_comp) * ((1 / v.eff_DCDC) - 1)) ...
                    - (P_HTS + x(2) + P_comp);

        P_H2_LHV_eff = x(1) / v.eff_FC;
        Q_dot_FC_stack = 0.95 * (P_H2_LHV_eff - x(1));
        P_TMS = 0.04 * Q_dot_FC_stack;

        eq2 = x(2) - P_TMS;

        F = [eq1; eq2];
    end

% Package outputs
P     = [P_M,P_DCAC,P_HTS, P_DCDC,P_FC_el,P_comp,P_TMS];
Q_dot = [Q_dot_M,Q_dot_DCAC,Q_dot_HTS,Q_dot_DCDC,Q_dot_FC_stack]; 

% Back out effective FC area from polarization (linearized Ucell = a*j + b)
apol= -0.2320;
bpol= 0.8956;
Ucell_avg= apol*jcell_max + bpol;
A_FC_eff= P_FC_el./(Ucell_avg*jcell_max); % [cm^2]

end

%% ========================================================================
%  PowerModel2j_IZEA_light
%  ------------------------------------------------------------------------
%  PURPOSE:
%    Similar to PowerModel2_* but takes A_FC_eff as input and returns the
%    mission-point current density jcell as output.
%
%  INPUTS:
%    Pprop, altitude, speed as before
%    A_FC_eff           [cm^2]
%    P_FC_el_initial    [W]
%    P_star_TMS         [W]
%
%  OUTPUTS:
%    P, Q_dot as before, plus:
%    m_dot_air [kg/s], jcell [A/cm^2]
%  ========================================================================
function [P,Q_dot,m_dot_air,jcell] = PowerModel2j_IZEA_light(Pprop,altitude,speed,A_FC_eff,P_FC_el_initial,P_star_TMS)

% Initial Parameters & efficiency maps
Pprop_initial = 16242000;

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

    v.eff_FC   = polyval(eff_fc_fit,   Pprop / Pprop_initial);
    v.eff_DCDC = polyval(eff_DCDC_fit, Pprop / Pprop_initial);
    v.eff_HTS  = polyval(eff_HTS_fit,  Pprop / Pprop_initial);
    v.eff_DCAC = polyval(eff_motorD_fit, Pprop / Pprop_initial);
    v.eff_M    = polyval(eff_motor_fit,  Pprop / Pprop_initial);
    
% First Calculation: propagate to HTS
components = simPar_IZEA_light(v,Pprop);
P_M=components(1);
Q_dot_M=components(2);
P_DCAC=components(3);
Q_dot_DCAC=components(4);
P_HTS=components(5);
Q_dot_HTS=components(6);

% Solve for FC & TMS
x0 = [P_FC_el_initial, P_star_TMS];
x_solution = fsolve(@target_fun, x0);
P_FC_el = x_solution(1);
P_TMS = x_solution(2);

% Compressor + m_dot using A_FC_eff → first compute jcell internally
[P_comp, m_dot_air] = PEMFCj(P_FC_el, altitude, speed, A_FC_eff);

% DC/DC, FC heat, TMS closure
Q_dot_DCDC = (P_HTS + (P_comp + P_TMS)) .*((1./v.eff_DCDC)-1);
P_DCDC = P_FC_el - Q_dot_DCDC;
P_H2_LHV_eff = P_FC_el / v.eff_FC;
Q_dot_FC_stack = 0.95 * (P_H2_LHV_eff - P_FC_el);
P_TMS = 0.04 * Q_dot_FC_stack;

    function F = target_fun(x)
        [P_comp, ~] = PEMFCj(x(1), altitude, speed, A_FC_eff);

        eq1 = x(1) - ((P_HTS + x(2) + P_comp) * ((1 / v.eff_DCDC) - 1)) ...
                    - (P_HTS + x(2) + P_comp);

        P_H2_LHV_eff = x(1) / v.eff_FC;
        Q_dot_FC_stack = 0.95 * (P_H2_LHV_eff - x(1));
        P_TMS = 0.04 * Q_dot_FC_stack;

        eq2 = x(2) - P_TMS;

        F = [eq1; eq2];
    end

P =  [P_M,P_DCAC,P_HTS,P_DCDC,P_FC_el,P_comp,P_TMS,P_H2_LHV_eff];
Q_dot = [Q_dot_M,Q_dot_DCAC,Q_dot_HTS,Q_dot_DCDC,Q_dot_FC_stack]; 

% Recover j from P_FC_el and A_FC_eff
jcell = jcalculation(P_FC_el,A_FC_eff);

end

%% ========================================================================
%  simPar_IZEA_light
%  ------------------------------------------------------------------------
%  PURPOSE:
%    Propagate power from propulsor upstream through DC/AC and HTS,
%    computing intermediate heats via efficiency models.
%
%  INPUTS:
%    v         - struct with fields eff_M, eff_DCAC, eff_HTS
%    Pprop [W] - demanded mechanical/electrical at propulsor node
%
%  OUTPUT:
%    components = [P_M, Q_M, P_DCAC, Q_DCAC, P_HTS, Q_HTS]
%  ========================================================================
function components=simPar_IZEA_light(v,Pprop)

P_Downstr = Pprop;

Q_dot_M = P_Downstr * ((1./v.eff_M)-1);
P_M = P_Downstr + Q_dot_M;

P_Downstr = P_M;

PF_M = 0.89;

Q_dot_DCAC = (P_Downstr./PF_M) * ((1./v.eff_DCAC)-1);
P_DCAC = P_Downstr + Q_dot_DCAC;

P_Downstr = P_DCAC;

Q_dot_HTS = P_Downstr * ((1./v.eff_HTS)-1);
P_HTS = P_Downstr + Q_dot_HTS;

components=[P_M  Q_dot_M  P_DCAC  Q_dot_DCAC  P_HTS  Q_dot_HTS];   

end

%% ========================================================================
%  PEMFC
%  ------------------------------------------------------------------------
%  PURPOSE:
%    Compressor power and air mass flow for a given FC electrical power and
%    target jcell (jcell_max). Uses standard compressible relations and an
%    ambient model vs altitude.
%
%  INPUTS:
%    P_FC_el [W], altitude [m], speed [m/s], jcell [A/cm^2]
%
%  OUTPUTS:
%    Pcomp         [W]   - compressor electrical power
%    m_dot_air_in  [kg/s]- air mass flow to FC cathode
%  ========================================================================
function [Pcomp,m_dot_air_in] = PEMFC(P_FC_el, altitude, speed, jcell)  % i have changed [Pcomp,Q_FCHX] by Pcomp
    % Constants
    ss = 343;  % Speed of sound [m/s]
    cp_air = 1000;  % [J/kg/K]
    y_air = 1.4;  % Gamma air
    Ma = speed / ss;
    M_air = 0.02896;      % [kg/mol]
    F = 96485;            % Faraday [C/mol]
    x_O2 = 0.21;          % [-]
    lamda_O2 = 1.8;       % stoich. excess air ratio
    eff_comp_m = 0.97;
    eff_comp_el = 0.94;
    eff_comp_pc = 0.95;
    eff_comp_s = 0.76;
    p_FCHX_in = 175000;   % [Pa]
    R_sp_air = 287;       % [J/kg/K]
    eff_pr = 0.75;
    % T_FCHX_out = 358; % K
    
   a=-0.2320;  % polarization slope [V·cm^2/A]
   b=0.8956;   % intercept [V]
   Ucell = a*jcell + b;      % average cell voltage [V]
   AFC = P_FC_el /(Ucell * jcell);  % required active area [cm^2]
    
    % Ambient conditions (ISA-like approximation)
    pamb = 101325 * (1 - 0.0065* (altitude./288.15))^5.2561;
    Tamb = 288.15-6.5*(altitude./1000)+24*(1-(altitude./12664)); % [K]

    % Mass flow calculation from stoichiometry
    m_dot_air_in = (jcell * AFC * lamda_O2 * M_air) / (4 * F * x_O2);
    
    % Compressor thermodynamics
    T_comp_in = Tamb + ((sqrt(y_air * R_sp_air * Tamb) * Ma)^2) / (2 * cp_air);
    p_comp_in_stat = pamb * (T_comp_in / Tamb)^(y_air / (y_air - 1));
    p_comp_in = eff_pr * (p_comp_in_stat - pamb) + pamb;
    delta_h_comp = (1 / eff_comp_s) * cp_air * T_comp_in * ((p_FCHX_in / p_comp_in)^(R_sp_air / cp_air) - 1);
    
    % Electrical power to drive compressor
    Pcomp = (delta_h_comp * m_dot_air_in) / (eff_comp_m * eff_comp_el * eff_comp_pc);
    % T_FCHX_in = T_comp_in + (delta_h_comp / 1000);
    % Q_dot_FCHX = cp_air * m_dot_air_in * (T_FCHX_in - T_FCHX_out);
end

%% ========================================================================
%  PEMFCj
%  ------------------------------------------------------------------------
%  PURPOSE:
%    Same as PEMFC, but jcell is computed from P_FC_el and A_FC_eff using
%    the polarization relation (via jcalculation).
%  ========================================================================
function [Pcomp,m_dot_air_in] = PEMFCj(P_FC_el, altitude, speed, AFC)  % i have changed [Pcomp,Q_FCHX] by Pcomp
    % Constants
    ss = 343;  % Speed of sound
    cp_air = 1000;  % J/kg*k
    y_air = 1.4;  % Gamma air
    Ma = speed / ss;
    M_air = 0.02896;
    F = 96485;
    x_O2 = 0.21;
    lamda_O2 = 1.8;
    eff_comp_m = 0.97;
    eff_comp_el = 0.94;
    eff_comp_pc = 0.95;
    eff_comp_s = 0.76;
    p_FCHX_in = 175000;
    R_sp_air = 287;
    eff_pr = 0.75;
    % T_FCHX_out = 358; % K
    
    % jcell from polarization (quadratic/iterative solver)
    jcell = jcalculation(P_FC_el,AFC);
    
    % Ambient conditions
    pamb = 101325 * (1 - 0.0065* (altitude./288.15))^5.2561;
    Tamb = 288.15-6.5*(altitude./1000)+24*(1-(altitude./12664)); % change
    
    % Mass flow calculation
    m_dot_air_in = (jcell * AFC * lamda_O2 * M_air) / (4 * F * x_O2);
    
    % Compressor calculations
    T_comp_in = Tamb + ((sqrt(y_air * R_sp_air * Tamb) * Ma)^2) / (2 * cp_air);
    p_comp_in_stat = pamb * (T_comp_in / Tamb)^(y_air / (y_air - 1));
    p_comp_in = eff_pr * (p_comp_in_stat - pamb) + pamb;
    delta_h_comp = (1 / eff_comp_s) * cp_air * T_comp_in * ((p_FCHX_in / p_comp_in)^(R_sp_air / cp_air) - 1);
    
    % Return values
    Pcomp = (delta_h_comp * m_dot_air_in) / (eff_comp_m * eff_comp_el * eff_comp_pc);
    % T_FCHX_in = T_comp_in + (delta_h_comp / 1000);
    % Q_dot_FCHX = cp_air * m_dot_air_in * (T_FCHX_in - T_FCHX_out);
end

%% ========================================================================
%  massCalc
%  ------------------------------------------------------------------------
%  PURPOSE:
%    Convert component rating requirements into mass estimates via simple
%    power-to-weight (or heat-to-mass) scaling laws, then build cumulative
%    stacked masses for plotting.
%
%  INPUTS: (units in W, kg/s, cm^2 as noted)
%    P_M, P_HTS, P_DCAC, P_DCDC, A_FC_eff, P_TMS, P_comp, m_dot_air
%
%  OUTPUTS:
%    m1..m9 cumulative masses [kg] in stacking order:
%    [stack, +comp, +hum, +BoP, +TMS, +DCDC, +HTS, +DCAC, +Motor]
%  ========================================================================
function [m1,m2,m3,m4,m5,m6,m7,m8,m9] = massCalc(P_M,P_HTS,P_DCAC,P_DCDC,A_FC_eff,P_TMS,P_comp,m_dot_air)%P_DCDC,P_FC_el,P_TMS

% n_SCM=3500;%rpm
% n_prop=1200;
% k_gear=0.21;
% m_gear=k_gear*((((P_SCM).^0.76)*((n_SCM).^0.13))/n_prop.^0.89);

%% Mass SCM (motor)
P_to_m_M =12000;
m_M= P_M * (1./P_to_m_M);

%% mass inverter constant
PTW_DCAC=27000; %W/kg
m_DCAC= P_DCAC/(PTW_DCAC);%10

%% mass HTS constant
P_to_m_HTS = 25000;
m_HTS=P_HTS * (1./P_to_m_HTS); % better estimate

%% mass converter
P_to_m_DCDC = 27000; %W/kg
m_DCDC = P_DCDC * (1./P_to_m_DCDC);

%% mass TMS
% Use optimal parameters as IZEA_light will use skin cooling
P_to_Q_TMS = 0.04;           % electrical power per unit heat removed
Q_dot_TMS= P_TMS./P_to_Q_TMS;
Q_to_m_TMS = 10000; % kj/s/kg   (scaling)
m_TMS = Q_dot_TMS./Q_to_m_TMS;
% m_TMS=P_TMS./10000;

%% mass stack,comp,hum,BoP non constant
rho_A_FC_eff=1.65;   % [kg per 10^4 cm^2] → m_stack = rho_A * A_FC_eff * 1e-4
% rho_comp=0.0004;
rho_hum=70;          % [kg / (kg/s)] for humidifier scaling
lamda_BoP=0.2;       % BoP fraction of (stack+comp+hum)

% m_stack = P_FC_el/2200;
m_stack= rho_A_FC_eff * A_FC_eff*1e-4;
% m_comp=rho_comp*P_comp;
m_comp=P_comp/11000;
m_hum=rho_hum * m_dot_air;
m_BoP= lamda_BoP*(m_stack+m_comp+m_hum);

%% Total mass return (cumulative stacking)
m1 = m_stack;
m2 = m1 + m_comp;
m3 = m2 + m_hum;
m4 = m3 + m_BoP;
m5 = m4 + m_TMS;
m6 = m5 + m_DCDC;
m7 = m6 + m_HTS;
m8 = m7 + m_DCAC;
m9 = m8 + m_M;

end

%% ========================================================================
%  jcalculation
%  ------------------------------------------------------------------------
%  PURPOSE:
%    Solve for j (A/cm^2) from P and A using linear polarization
%    Ucell = a*j + b, resulting quadratic in j. Uses closed-form solution
%    with fallback to bisection if needed.
%  ========================================================================
function jcell = jcalculation(P, A)
    %% Parameters
    a = -0.2320; % V·cm^2/A
    b = 0.8956;  % V
    
    %% Calculation using more efficient root-finding
    % Solve: A*((a*j^2 + b*j)/P) = 1  → a*j^2 + b*j - P/A = 0
    c = -P/A;
    discriminant = b^2 - 4*a*c;
    
    if discriminant < 0
        % No real roots
        jcell = NaN;
        return;
    end
    
    % Quadratic formula
    j1 = (-b + sqrt(discriminant))/(2*a);
    j2 = (-b - sqrt(discriminant))/(2*a);
    
    % Choose the positive, feasible root (≤ 1.6 A/cm^2)
    if j1 > 0 && j1 <= 1.6
        jcell = j1;
    elseif j2 > 0 && j2 <= 1.6
        jcell = j2;
    else
        jcell = NaN;
    end
    
    % Quick consistency check
    sum = A*((a*jcell^2 + b*jcell)/P);
    error = abs(sum - 1);
    
    % Fallback to iterative method if needed
    if error > 1e-2 || isnan(jcell)
        jcell = iterative_solution(P, A, a, b);
    end
end

%% ========================================================================
%  iterative_solution (helper)
%  ------------------------------------------------------------------------
%  PURPOSE:
%    Bisection fallback for j calculation if quadratic root is unsuitable.
%  ========================================================================
function j = iterative_solution(P, A, a, b)
    % Initial guess
    j = 1; % mid-range
    
    % Bisection bounds and tolerance
    j_min = 0;
    j_max = 1.6;
    error_threshold = 1e-4;
    
    for i = 1:100 % iteration cap
        sum = A*((a*j^2 + b*j)/P);
        error = sum - 1;
        
        if abs(error) < error_threshold
            break;
        end
        
        % Update bounds
        if error > 0
            j_max = j;
        else
            j_min = j;
        end
        
        % Bisection
        j = (j_min + j_max)/2;
        
        % Safety break
        if j_max - j_min < 1e-6
            break;
        end
    end
end
