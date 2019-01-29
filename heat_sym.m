clear all
clc

%how much water do we need if the outside is at it's coldest for 7 days?
%how much heat would be needed if there was a pipe circulating the house?
%assume insulation amount

V_side1 = 3; %[m]
V_side2 = 3; %[m]
V_floor = 3; %[m]
V_house = V_side1 * V_side2 * V_floor; %[m^3 = 10ft^3]

%Insulation properties 
%http://highperformanceinsulation.eu/wp-content/uploads/2016/08/Thermal_insulation_materials_made_of_rigid_polyurethane_foam.pdf
k_ins_wall = 0.025; %Rigid polyurethane foam (PUR/PIR) insulation; [(W/(mK)])
rho_ins = 30; %density [kg/m^3]
ksp_ins = 1.5; %specific heat [kJ/(kgK)]
A_ins = V_side1 * V_side2; %insulation Area [m^2]
t_ins_wall = .30; %thickness; [m]

%water drum; https://www.uline.com/Product/Detail/S-14364/Drums/Steel-Drum-with-Lid-30-Gallon-Open-Top-Unlined?pricode=WB0299&gadtype=pla&id=S-14364&gclid=Cj0KCQiA7briBRD7ARIsABhX8aBXzahTKj4O7YuUj8KMNm21OF2ZJds-dQyWbdpo2vGf0dmV7Cs9M0kaAoHpEALw_wcB&gclsrc=aw.ds
%drum temperature range: 233.15K to 435.928K (-40F to 325F)
V_drum = 30 * 0.00378541; % [gal] * [m^3/gal] = [m^3]
t_drum = 1.1E-3; %[m]
k_steel = 385.0; %[W/mK]
k_ins_drum = 0.04; %[W/mK] fiberglass http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/thrcn.html
r_drum_out = 19 * 2.54 / 100 / 2; %[in] * [cm/in] * [m/cm] = [m]
r_drum_in = 18 * 2.54 / 100 / 2; %[in] * [cm/in] * [m/cm] = [m]
r_drum_ins = (r_drum_out + 2) * 2.54 / 100 / 2; %[in] * [cm/in] * [m/cm] = [m]
h_drum = 28.75 * 2.54 / 100; %[in] * [cm/in] * [m/cm] = [m]
SA_drum_in = pi * 2 * r_drum_out * h_drum;
t_ins_drum = .025; %[m]

T_inside = 55 + 273; %[K]
T_outside = -36.7 + 273; %[K]

%find rate that heat leaves greenhouse
Qdot_out = 6 * k_ins_wall * A_ins * (T_inside - T_outside) / t_ins_wall; %Heat transfered per second; [W]

T_water = 366.483; %[K]

%assume: rate that heat is being replaced = rate of heat leaving the
%greenhouse (steady-state) AKA Qdot_out = Qdot_in
Qdot_in = Qdot_out;

%Heat radiating off insulated water drum
R_wall = log(r_drum_out / r_drum_in) / (2 * pi * k_steel * h_drum);
R_insulation = log(r_drum_ins / r_drum_out) / (2 * pi * k_ins_drum * h_drum);
R_drum_tot = R_wall + R_insulation;
Q_drum_rad_out = (T_water - T_inside) / R_drum_tot;

t_storage = 7 * 24 * 3600; %[days] * [24hrs/1day] * [3600s/1hr] = [s]

% find required energy to keep room at same temp
% lose some heat from radiation, but can add that to heating the room

Q_storage = Qdot_in * t_storage; %[J/s] * [s] = [J]





