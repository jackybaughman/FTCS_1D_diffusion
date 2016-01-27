%heat1Dexplicit.m
%
% Solves the 1D heat equation with an explicit finite difference scheme
% FTCS
% Ts abrupt change
% JSB w/ help from USC code from GEOL557 Jan. 2016

clear all
figure(1)
clf

%% Physical parameters
L = 400; % length of modeled domain [m]
Ts = -6; % Surface T [C]
kappa = 1e-6; % Thermal diffusivity of rock [m2/s]
k = 2; % conductivity
Q = .035;% heat flux
rho = 2000;
Cp = 2000;
day = 3600*24; % # seconds per day
dt = 10*day; % Timestep [s]


%% Numerical parameters
nx = 201; % Number of gridpoints in x-direction
nt = 200; % Number of timesteps to compute
dz = L/(nx-1); % Spacing of grid
z = 0:dz:L;% Grid
% Setup initial temperature profile
T = ones(size(z))*Ts + ((Q*z)/k);
time = 0;

%% 
%%Part B
% analytic solutions animated

Tsurface = -3;
Tnought = -6;
Tideal = T;
Tideal = T +(Tsurface-Tnought)*((1-erf(z./(2*sqrt(kappa * dt*nt)))));
figure(1)
plot(T,z,'r','linewidth',2)
  xlabel('Temperature (C)','fontname','arial','fontsize',21)
    ylabel('Depth (m)','fontname','arial','fontsize',21)
    set(gca,'fontsize',18,'fontname','arial')
    set(gca,'YDIR','reverse')
    hold on
    
%% Part C: File:  capethompson.m
%
%  Load data from capethompson.dat and plot it with symbols

load capethompson.dat;         %  read data into PDXprecip matrix 
depth = capethompson(:,1);     %  copy first column of PDXprecip into month
Temperature = capethompson(:,2);    %  and second column into precip


%% run
for n=1:nt % Timestep loop
% Compute new temperature
Tnew = zeros(1,nx);
for i=2:nx-1
Tnew(i) = T(i)+(((T((i+1))-2*T(i)+T((i-1)))/(dz*dz))*(dt*kappa));
end
% Set boundary conditions
Tnew(1) = -3;
Tnew(nx) = T(nx);
% Update temperature and time
T = Tnew;
time = time+dt;
% Plot solution
figure(1), clf
plot(Tnew,z);
xlabel('Temperature [^oC]')
ylabel('x [m]')
set(gca,'YDIR','reverse')
set(gca,'fontsize',14,'fontname','arial')
title(['Temperature evolution after ',num2str(time/day),' days'])
hold on
plot(Tideal,z,'r','linewidth',2)
hold on 
plot(Temperature, depth,'o');     %  plot temp vs. depth with circles
ylabel('depth (m)');              %  add axis labels and plot title
xlabel('Temperature (C)');
set(gca,'YDIR','reverse')
hold off

drawnow
end


