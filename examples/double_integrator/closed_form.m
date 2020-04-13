clc
clear
dt =0.1;
time = 0:dt:10;
x0 = [1; 1];
tf = 10;
u = -2/tf^2 * ( 3*x0(1) +2*x0(2)*tf) + 6/tf^3 * (2*x0(1) + x0(2)*tf)* time;

x = ones(2,length(time));
for i=1:length(time)-1
    k1=state_dynamics(x(:,i),u(i));
    k2=state_dynamics(x(:,i)+0.5*k1*dt,u(i));
    k3=state_dynamics(x(:,i)+0.5*k2*dt,u(i));
    k4=state_dynamics(x(:,i)+dt*k3,u(i));
    x(:,i+1) = x(:,i)+1/6*dt*(k1+2*k2+2*k3+k4);
end
%[~,x] = ode45(@(t,x) dynamics(t,x,time,u), time, x0);

plot(time,u)
figure
plot(time,x(1,:))
%plot(time,x(:,1))
function dx= state_dynamics(x,u)
    dx(1,1) = x(2);
    dx(2,1) = u;
end
function dx = dynamics(t,x,time,u)
    dx(1,1) = x(2);
    dx(2,1) = interp1(time,u,t);
end