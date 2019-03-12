function [waterTempOut] = deltaWater(waterTempIn)
T1 = 95+273;
T2 = 50+273;
slope = 7/(T1-T2);
b = 12-slope*T1;
waterTempOut = slope*waterTempIn + b;
end

