clear all
clc
iters = 200;
repeat = 30;
maxx = zeros(iters,1);
avgmax = zeros(repeat,1);

Kb =1.38e-23;   %Bolztman const
T =300000;         %Temperature
R = 3.5e-9 ;    %Radius
m = 1.54e-11;   %Mass
mu = 0.001;     %Viscosity
lambda = 6*pi*mu*R; %Stokes law
D = lambda*Kb*T ; %Einstein relation

%Initial conditions and parameters
t = 1;   
j= 200;
N = 10^4;    
dt = t/N;
for k = 1:repeat
    
for z=1:iters
p = zeros(j,N);
x = zeros(j,N);
p(:,1) = 0;

%Calculate Brownian forces
dW = sqrt(dt*2*D)*randn(j,N);

%Integrate
for i=2:N
        p(:,i) = p(:,i-1) - dt*(lambda)*p(:,i-1)/m + dW(:,i);
        x(:,i) = x(:,i-1) + (p(:,i)/m)*dt;
end;

maxx(z,1) = max(x(:,N));
avgmax(k,1) = mean(maxx);
end;
end;
plot(1:repeat,avgmax)
