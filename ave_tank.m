%inputs

clear all


%water drum; https://www.uline.com/Product/Detail/S-14364/Drums/Steel-Drum-with-Lid-30-Gallon-Open-Top-Unlined?pricode=WB0299&gadtype=pla&id=S-14364&gclid=Cj0KCQiA7briBRD7ARIsABhX8aBXzahTKj4O7YuUj8KMNm21OF2ZJds-dQyWbdpo2vGf0dmV7Cs9M0kaAoHpEALw_wcB&gclsrc=aw.ds
%drum temperature range: 233.15K to 435.928K (-40F to 325F)
t_d = 1.1E-3; %[m]
k_d = 385.0; %[W/mK]
k_ins = 0.04; %[W/mK] fiberglass http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/thrcn.html
t_ins = .012; %[m]
ar_d = .633; %aspect ratio drum

T_outside = -38.16 + 273; %[K]

%Water properties
T_w = 35 + 274; %[K]
rho_w = 963.33; %kg/m^3 https://water.usgs.gov/edu/density.html
c_w = 3.8204*1000; % [kJ/kgK]*[J/kJ] = [J/kgK] Isochoric heat capacity of water at 363K https://www.engineeringtoolbox.com/specific-heat-capacity-water-d_660.html
d_T = T_w - T_outside;
y = 7 * 24 * 3600; %seconds in a week

gal = 1;
t_req = 7 * 24 * 3600; % 7 days [s]
t_storage = 0; %initial t_storage

T_w_0 = 35 + 273; %[K]
T_gh = 15.5 + 273; %[K]

u_w = 2; %[m/s]
r_pipe = .005; %[m]
A_pipe = r_pipe^2 * pi;
m_dot = A_pipe * rho_w * u_w; %mass flow rate [kg/s]

t_step = 1; %[s]
m_in = m_dot * t_step; %[kg] kg water in for a step

%while loop running until the maximum time that water can be stored > time required by the GH
%step 1 - amount of gallons
%step 2 - find new mass and volume of water
%step 3 - for as long as the temperature of the water is greater than the
%greenhouse...
%step 3a - find new water temp given heat loss through tank
%step 3b - find new average water temp given a unit of cold in and hot out
%step 3c - find the time it takes to cool down

while t_storage < t_req 
    
    %find tank volume and associated variables
    gal = 1 + gal;
    conv = 0.00378541; %1 gal = 0.00378541 m^3
    V_w = gal * conv; %[m^3]; 1 gal = 0.00378541 m^3
    h_d = (4 * V_w / (pi * ar_d^2))^(1/3); % heigth drum [m]
    r_d = ar_d * h_d / 2; %radius drum [m]
    m_w = V_w * rho_w; %[kg]
    
    %tank cooling down
    T_w = T_w_0;
    
    while T_w > T_gh
        
        t_storage = t_storage + 1; %seconds elapsed
        
        %find T_w after losses in the period of 1 second
        R_d = log((r_d + t_d)/r_d) / (2 * pi * k_d * h_d); %thermal resistance through drum
        R_ins = log((r_d + t_d + t_ins)/(r_d + t_d)) / (2 * pi * k_ins * h_d); %thermal resistance through drum insulation
        Qdot_loss = (T_w - T_outside) / (R_d + R_ins); %[W]
        Q_loss = Qdot_loss * t_step; %energy loss for the time period [J]
        Q_storage = m_w * c_w * (T_w - T_outside); %energy in tank [J]
        Q_tot = Q_storage - Q_loss; %energy in tank after loss [J]
        T_w = Q_tot / (m_w * c_w) + T_outside; %new T_w after energy loss through tank wall [K]
        
        hot_count = m_w - m_in; %[kg]
        T_w = (hot_count * T_w + m_in * T_gh) / m_w; %find new average tank temperature [K]
        
    end
    
end

gal %how many gallons are required to last the time period required?
r_d
h_d
V_w






