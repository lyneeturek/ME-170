function [q_loss] = qLossGH(T_air)
%[INPUT] Dimensions of structure (ft)
dims=[20,10,10];

%[INPUT] Thermal values
T_inf=-20; %Outside air temperature (C)
T_air=18; %Optimal cabbage-growing temperature (C)
th=6.5; %Wall thickness (in)
    %SIP with R-value of 24
K=1/25.6005; %Conductivity of insulation (W/mK)
windspeed=4; %Wind speed (m/s) - cut in wind speed for turbines
F_s=0.7145; %Shape factor for conduction

%Constant calculations
dims=dims*unitsratio('meters','feet');
th_metric=th*unitsratio('meters','in');
A=2*dims(1)*dims(2)+2*dims(1)*dims(3)+2*dims(2)*dims(3);
roughness=[8.23,4.0,-0.057]; %Surface roughness coefficients (clear pine)
h=roughness(1)+roughness(2)*windspeed+roughness(3)*windspeed^2;

%Thermal resistances
R_cond=th_metric/(K.*A)*F_s;
R_conv=1/(h.*A);

q_loss=(T_air-T_inf)./(R_cond+R_conv);
end

