m = 1;
n = 2;
p = 2;
X0 = [ 1, 1];
Yf = [0; 0];
R = [1];
t0 = 0;
tf = 10;
dt = 0.1;
tolerance = 0.02;
N = floor(abs(tf - t0)/dt) + 1;
U_guess = ones(m,N);
maxIter = 200;
%[t, u, x] = GMPSP(@forward,@output, @dF_dX, @dF_dU, @dh_dX, 0 , 10,...
%    0.1, X0, Yf, R, ones(1,100),0.02,2000,1);
[t, X_opt, U_opt] = GMPSP(@forward,@output, @dF_dX, @dF_dU, @dH_dX, t0 , tf,...
    dt, X0, Yf, R, U_guess, tolerance, maxIter);

function dx = forward(t,x,u)
dx = [x(2); u];
end

function y = output(x)
y = [x(1) ; x(2)];
end

function out = dF_dX(x,u)
out = [0 1 ; 
       0 0];
end

function out = dF_dU(x,u)
out = [0 ; 1];
end

function out = dH_dX(x,u) 
out = [1 0; 0 1];
end