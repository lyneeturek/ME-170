close all; clear all; clc;
% Product Info
% Vinyl Hose -> https://www.webstaurantstore.com/manitowoc-ice-4420383-tubing-clear-vinyl-5-8id-7-8od/HP4420383.html?utm_source=Google&utm_medium=cpc&utm_campaign=GoogleShopping&gclid=EAIaIQobChMI7aafx8DI4AIVhCCtBh1yugwjEAQYASABEgIUT_D_BwE
% Clay Pot -> https://www.homedepot.com/p/Pennington-12-in-Terra-Cotta-Clay-Pot-100043019/100333337
k_ins = 0.033;%[W/mK] https://www.engineeringtoolbox.com/thermal-conductivity-d_429.html
k_vinyl = 0.25;%[W/mK] https://www.engineeringtoolbox.com/thermal-conductivity-d_429.html
ins_thickness = 6 * 0.0254;%[in]*[m/in] = [m]
U_air = 22;%[m/s] https://www.reviews.com/hair-dryer/ (used lowest airspeed of high-end blow dryer and subtracted 10 mph to get that of our budget blow dryer)

pebbles_percent = 0.7;
air_percent = 1-pebbles_percent;
m_pebbles = 62000;%kg
rho_pebbles = 2650;%[kg/m^3] used sandstone
rho_air = 1.225;%[kg/m^3] air density at sea level
rho_pebble_air = pebbles_percent*rho_pebbles + air_percent*rho_air;

cp_pebbles = 920;%[J/kgK] Specific heat of sandstone https://www.engineeringtoolbox.com/specific-heat-solids-d_154.html
L = 2.9;%[m]

t_final = 60*60*24*4;%4 days in seconds [s]

T_inlet = (17 + 273.15);%[K]
T_pebbles_initial = (95 + 273.15);%[K]
T_ambient = -20 + 273.15;%[K]

inner_hose_radius = 2 * 0.0254;%[meters]
outer_hose_radius = 3 * 0.0254;%[meters]
inner_pebbles_radius = outer_hose_radius;
outer_pebbles_radius = (((m_pebbles/pebbles_percent)/rho_pebble_air*pi*L)+outer_hose_radius^2)^(1/2);
inner_ins_radius = outer_pebbles_radius;%[meters]
outer_ins_radius = inner_ins_radius + ins_thickness;%[meters]
pipe_surface_area = 2 * pi * inner_hose_radius * L;%[m^2] area of hose for convection
A_cross_section = pi * inner_hose_radius^2;%[m^2]

cp_air = 1005;%[J/KgK]
h_air = 10.45 - U_air + 10 * U_air^(1/2);% for 2 < U_air < 20 https://www.engineeringtoolbox.com/convective-heat-transfer-d_430.html
Nu_tube = 3.66;% assuming fully developed laminar flow in tube http://www.me.nchu.edu.tw/lab/lab516/2014/Heat%20Transfer-PDF-Incropera-1/8b.pdf

R_hose = log(outer_hose_radius/inner_hose_radius)/(2*pi*L*k_vinyl);%[K/Watt]
R_ins = log(outer_ins_radius/inner_ins_radius)/(2*pi*L*k_ins);%[K/Watt]
R_convection = 1/(h_air*pipe_surface_area);%[K/Watt]

n = 10000;% number of points
t = linspace(0,t_final, n);
dt = t(2)-t(1);

%assume y distribution of temperature in pebbles is constant
m_dot_air = rho_air*U_air*A_cross_section;

%initialize vectors and evaluate initial conditions
T_pebbles = zeros(1, length(t));
T_pebbles(1) = T_pebbles_initial;
T_air = zeros(1, length(t));
T_air(1) = T_inlet;
Q_dot_loss_pebbles = zeros(1, length(t));
Q_dot_out_pebbles = zeros(1, length(t));
Q_dot_loss_pebbles(1) = (T_pebbles_initial-T_ambient)/(R_ins);
Q_dot_out_pebbles(1) = (T_pebbles_initial-T_inlet)/(R_hose + R_convection);


for i = 2:length(t)
Q_dot_out_pebbles(i) = (T_pebbles(i-1) - T_air(i-1))/(R_hose + R_convection);
Q_dot_loss_pebbles(i) = (T_pebbles(i-1) - T_ambient)/R_ins;

Q_dot_in_air = Q_dot_out_pebbles(i);

cp_air = 1.9327e-10 * T_air(i-1)^4 - 7.9999e-7 * T_air(i)^3 + 1.1407e-3 * T_air(i)^2 - 4.4890e-1 * T_air(i) + 1.0575e3;% [J/kgK] file:///Users/alejandroballesteros/Downloads/air_cp_plot.pdf
dT_air = Q_dot_in_air/(m_dot_air*cp_air);
dT_pebbles = (Q_dot_loss_pebbles(i) + Q_dot_out_pebbles(i))*dt/(m_pebbles*cp_pebbles);

T_pebbles(i) = T_pebbles(i-1) - dT_pebbles;
T_air(i) = T_air(i-1)+ dT_air;
end

figure(1)
plot(t/(3600*24), T_pebbles-273.15)
hold on
plot(t/(3600*24), T_air-273.15);
title('Air and Pebble Temperature (C) vs Time (min)')
xlabel('Time (days)');
ylabel('Temperature (C)');
legend('Temperature of Pebbles', 'Temperature of Air');


figure(2)
plot(t/(3600*24), Q_dot_out_pebbles);
hold on
plot(t/(3600*24), Q_dot_loss_pebbles);
title('Rate of Heat Transfer (Watts) vs Time (min)');
xlabel('Time (days)');
ylabel('Qdot (Watts)');
legend('Heat Transfer from Pebbles to Air', 'Heat Transfer from Pebbles to Surroundings');
