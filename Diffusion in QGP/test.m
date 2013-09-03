clear all
clc


Kb =1.38e-23;   %Bolztman const
T =300000;         %Temperature
R = 3.5e-9 ;    %Radius
m = 1.54e-11;   %Mass
mu = 0.001;     %Viscosity
lambda = 6*pi*mu*R; %Stokes law
D = lambda*Kb*T ; %Einstein relation

%Initial conditions and parameters
t = 2;   
j= 500;
N = 10^4;    
dt = t/N;

p = zeros(j,N);
x = zeros(j,N);
p(:,1) = 4e-14;

%Calculate Brownian forces
dW = sqrt(dt*2*D)*randn(j,N);

%Integrate

    for i=2:N
        p(:,i) = p(:,i-1) - dt*(lambda)*p(:,i-1)/m + dW(:,i);
        x(:,i) = x(:,i-1) + (p(:,i)/m)*dt;
    end;
   
Mp = mean(p);
Mp2 = mean((p/m).^2);
Mx = mean(x);
taxis = dt*(1:N);
figure()
subplot(211)
plot(taxis,Mp)
hold on
plot(taxis,p(1:3,:),'--r')
plot(taxis,p(1,1)*exp(-(lambda/m)*(dt*(1:N))),'r')
subplot(212)
plot(taxis,Mx)
hold on
plot(taxis,x(1:5,:),'--r')
%figure()
%plot(taxis,p)
%plot(taxis,x)
figure()
plot(taxis,Mp2)
hold on
plot(taxis,Kb*T/m)