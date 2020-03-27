X0 = [ 1, 1];
Yf = [0; 0];
R = [1];
[t, u, x] = MPSP(@forward,@output, @dF_dX, @dF_dU, @dh_dX, 0 , 10,...
    0.1, X0, Yf, R, ones(1,100),0.02,2000,1);

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

function out = dh_dX(x,u) 
out = [1 0; 0 1];
end