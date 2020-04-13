clc 
clear
syms r phi lam V GM gamma  
lam = r * V^2 / GM;
output = r/ ( (1- cos(phi))/(lam * cos(gamma)^2) + (cos(phi+gamma))/(cos(gamma)));

simplify(diff(output,r))
simplify(diff(output,V))
simplify(diff(output,gamma))
simplify(diff(output,phi))