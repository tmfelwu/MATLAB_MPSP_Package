clc
clear all
A=[0,1;0,0];
B=[0;1];
x_init=[1;1];
x_f=[0;0];
t=0:0.01:10;
u=ones(1,length(t));
x=zeros(2,length(t));
x(:,1)=x_init;
w11=ones(1,length(t));
w12=zeros(1,length(t));
w21=zeros(1,length(t));
w22=ones(1,length(t));

for i=1:length(t)
        w12(i)=-(t(i)-t(end));
end

w=[w11, w12; w21,w22];
bs=[w12;w22];
A_lambda=[-(-10)^3/3, 100/2;100/2, 10];
temp_b=zeros(2,length(t));
for i=1:length(t)-1
    k1=0.01*state_dynamics(x(:,i),u(i));
    k2=0.01*state_dynamics(x(:,i)+0.5*k1,u(i));
    k3=0.01*state_dynamics(x(:,i)+0.5*k2,u(i));
    k4=0.01*state_dynamics(x(:,i)+k3,u(i));
    x(:,i+1) = x(:,i)+1/6*(k1+2*k2+2*k3+k4);
end
iter=0;
while(x(1,end)^2+x(2,end)^2>1e-3)
    iter=iter+1;
    for i=1:length(t)
        temp_b(:,i)=[ -(t(i)-t(end))*u(i);u(i)];
    end
    b_lambda=trapz(t,temp_b,2);
    for m=1:length(t)
        u(m)=-transpose(bs(:,m))*((A_lambda)\(x(:,end)-b_lambda));
    end
    for i=1:length(t)-1
    k1=0.01*state_dynamics(x(:,i),u(i));
    k2=0.01*state_dynamics(x(:,i)+0.5*k1,u(i));
    k3=0.01*state_dynamics(x(:,i)+0.5*k2,u(i));
    k4=0.01*state_dynamics(x(:,i)+k3,u(i));
    x(:,i+1) = x(:,i)+1/6*(k1+2*k2+2*k3+k4);
    end
end
plot(t,x(1,:))
    