syms r V gamma phi delta mass thrust g
f = [ V* sin(gamma);
    thrust/mass * cos(delta) - g * sin(gamma);
    -thrust/(mass*V) * sin(delta) + (V/r - g/V)* cos(gamma);
    - V*cos(gamma)/ r];

% calculate df_dx
simplify(diff(f,r))
simplify(diff(f,V))
simplify(diff(f,gamma))
simplify(diff(f,phi))

% calculate df_du
simplify(diff(f,delta))