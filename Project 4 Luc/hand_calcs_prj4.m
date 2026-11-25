%Drone arm calculation verification
clear all; close all; clc;

%units/constants
intom = 0.0254;
lbftokg = 0.45359237;
g = 9.80665;


%geometry controls
do = 0.04; %1.5*intom; %in -> m, outer diameter
t = 0.0025; %0.0625*intom; %i-> m, wall thickness
l = 0.6096; %24*intom; %in -> m, arm length

%mech properties controls
E =  64*10^9; %Pa
v = 0.33; %poissons
rho = 1520; %kg/m3 density AL 6061 
G = E/(2*(1*v)); %shear mod

%loading controls
T = 30*lbftokg; %lbf -> kg
motor_mass = 1.5*lbftokg; %lbf -> kg

%calc geometry
area = (pi/4)*(do^2-(do-2*t)^2); %c/s area
I = (do^4-(do-2*t)^4)*(pi/64); %area MOI
J = (do^4-(do-2*t)^4)*(pi/32); %polar MOI
tube_weight = area*rho*l*g; %N

%arm bending and torsion
tipload = (T-motor_mass)*g;
bendingMoment = tipload*l - tube_weight*(l/2);
bendingStress = bendingMoment*(do/2)/I;
torsionalStress = T*(do/2)/J;

%arm deflection (max at tip)
bendingDeflection = tipload*l^3/(3*E*I) - tube_weight*l^3/(8*E*I);


%frequency 
mu = rho*area;
eta = motor_mass/(mu*l);
characteristic = @(b) 1 + cos(b)*cosh(b) + eta*b*(cos(b)*sinh(b) - sin(b)*cosh(b)); %uniform cantilever with point mass at tip vibrational characteristic equation - NACA TN 3459 p 17, eq13 and 
                                                                                    %https://userweb.eng.gla.ac.uk/sondipon.adhikari/fulltext/journal/ft51.pdf#page=15
fn = zeros(1,5);
for n = 1:5
    b = fzero(characteristic, [(n-1)*pi, n*pi]);
    fn(n) = b^2/(2*pi*l^2) * sqrt(E*I/mu);
end

q_V2 = -tube_weight/l; 

V2_midspanDeflection = 5*q_V2*l^4/(384*E*I);

%report findings
fprintf('Tipload: %.10f N\nCross Sectional Area: %.10e m2\nArea MOI: %.10e m4\nPolar MOI: %.10e m4\n',tipload,area, I, J)
fprintf('Tube weight: %.10f N\nBending stress: %.10f MPa\n', ...
    tube_weight, bendingStress/1e6);
fprintf('Bending deflection: %.10f mm\n',bendingDeflection*1e3);
fprintf('Natural Frequency Mode 1: %.10f Hz\nNatural Frequency Mode 2: %.10f Hz\nNatural Frequency Mode 3: %.10f Hz\nNatural Frequency Mode 4: %.10f Hz\nNatural Frequency Mode 5: %.10f Hz\n', fn(1),fn(2),fn(3),fn(4),fn(5))

fprintf('Midspan deflection: %.13f mm\n', V2_midspanDeflection*1e3);
