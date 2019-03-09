clear all
close all
clc

%% Parameters

%water drum; https://www.uline.com/Product/Detail/S-14364/Drums/Steel-Drum-with-Lid-30-Gallon-Open-Top-Unlined?pricode=WB0299&gadtype=pla&id=S-14364&gclid=Cj0KCQiA7briBRD7ARIsABhX8aBXzahTKj4O7YuUj8KMNm21OF2ZJds-dQyWbdpo2vGf0dmV7Cs9M0kaAoHpEALw_wcB&gclsrc=aw.ds
%allowable drum temperature range: 233.15K to 435.928K (-40F to 325F)
t_d = 1.1E-3; %[m]
k_d = 385.0; %[W/mK]
k_ins = 0.04; %[W/mK] fiberglass http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/thrcn.html
t_ins = 4 * .0254; %.012; %[m]
h_d = 9.5 * 12 * 2.54 / 100;%[m] = [ft] [12in/ft] [2.54cm/in] [1m/100cm]
n_nodes = 10;

%Water properties
rho_w = 963.33; %kg/m^3 https://water.usgs.gov/edu/density.html
c_w = 4.19*1000; % [kJ/kgK]*[J/kJ] = [J/kgK] Isochoric heat capacity of water at 363K https://www.engineeringtoolbox.com/specific-heat-capacity-water-d_660.html
k_w = .606;

%Temperatures
T_outside = -20 + 273; %[K] coldest temperature Kong could reach
T_w_0 = 95 + 273; %[K] intial water temperature
T_w_min = 50 + 273; %[K]
T_gh = 18 + 273; %[K] greenhouse temperature
T_gh_min = 10 + 273; %[K]
T_gh_max = 26 + 273; %[K]

%Greenhouse
conv = 0.00378541; %1 gal = 0.00378541 m^3
Qdotloss_gh = 1153; %[W]
Qloss_gh = Qdotloss_gh * 3600; %[J/hr]

%Vector Outputs
gal = linspace(1,200);
T_old = zeros(1,n_nodes);
[i1, i2] = size(gal);

%% Body

for j = 1:i2
    
    %reset counters and temperatures
    T_w = T_w_0;
    t_s = 0; %storage time
    T_old(1:n_nodes) = T_w; %old node temperature
    T = T_old;

    %new amount of water (step up)
    V_w = gal(j) * conv; %[m^3]; 
    r_d = sqrt(V_w / (pi * h_d)); %radius drum [m]
    m_w = V_w * rho_w; %[kg] mass water
    m_n = m_w/n_nodes; %[kg] mass of a node
    h_n = h_d/n_nodes; %[m] height of a node
    
    while T(n_nodes) > T_w_min
        
        t_s = t_s + 1; %time step up one
        
        %find the new mass of water required and period for the duty cycle
        %given the new outlet tank temperature
        dT_w = find_deltaT_w(T(n_nodes)); %find the water temperature drop given the water outlet temperature
        m_w_req = Qloss_gh / (c_w * dT_w); % kg required for an hour
        Vdot = 3 * 60 * conv; %[m^3/hr] = [gal/min] [60min/1hr] [.00378541 m^3/gal]
        mdot = Vdot * rho_w; %[kg/hr] = [m^3/hr] [kg/m^3]
        
        %how long does it take gh to cool down 3C
        t_ghcool % = some number;
        
        %will get x amount of energy out for every loop
        %keep running the pump until 
        
     
    %first node 
    T(1)  = ((m_n - mdot) * T_old(1) + mdot * T_gh) / m_n;

    %mass displacement into other nodes
        for i = 2:n_nodes
            T(i) = ((m_n - m_in) * T_old(i) + m_in * T_old(i-1))/ m_n;
        end

    T_old = T;

    
    
    
    
    
    end
    
    
    
    
    
    
    
    
    
    
    
    
end






