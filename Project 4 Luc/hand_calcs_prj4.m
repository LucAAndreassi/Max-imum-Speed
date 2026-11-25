%Drone arm calculation verification
clear all; close all; clc;

%units/constants
intom = 0.0254;
lbftokg = 0.4536;
g = 9.80665;


%geometry controls
do = 0.055285; %1.5*intom; %in -> m, outer diameter
t = 0.00201; %0.0625*intom; %i-> m, wall thickness
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


%report findings
fprintf('Tipload: %.2f N\nCross Sectional Area: %.6e m2\nArea MOI: %.6e m4\nPolar MOI: %.6e m4\n',tipload,area, I, J)
fprintf('Tube weight: %.2f N\nBending stress: %.2f MPa\n', ...
    tube_weight, bendingStress/1e6);
fprintf('Bending deflection: %.3f mm\n',bendingDeflection*1e3);
fprintf('Natural Frequency Mode 1: %.2f Hz\nNatural Frequency Mode 2: %.2f Hz\nNatural Frequency Mode 3: %.2f Hz\nNatural Frequency Mode 4: %.2f Hz\nNatural Frequency Mode 5: %.2f Hz\n', fn(1),fn(2),fn(3),fn(4),fn(5))


