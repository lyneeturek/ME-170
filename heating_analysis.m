clear all; close all; clc;

Qdot = 98000;%[J/s] assuming 98kwh every hour
dT = 45;%K (from 50C to 95C)
cp_pebbles = 945.5;%J/KgK
m_pebbles = 62000;%Kg
Q = m_pebbles*cp_pebbles*dT;
t = Q/Qdot;
t = t/(3600)