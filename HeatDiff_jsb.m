% Geothermal gradient - relaxing after major erosion event
% 
% FTCS
% 
% JSB code update from 
clear all
figure(1)
clf

%% Constants
Tsnaught = 20; %initial surface T C
Tstep = 5; % Temperature step change

k = 2; % conductivity
rho = 2000; % density, kg/m^3
c = 2000; % Specific heat capacity, J/kg*K
kappa = k/(rho*c); % Thermal diffusivity, m^2/s

Qm = .04; % mantle heat flux

%% Arrays
% z array
dz = 1; % z step m
zmax = 1000; % zmax m
z = 0:dz:zmax;% z array

% x array
dx = 1; % x step m
xmax = 1000; % xmax m
x = 0:dz:xmax;%x array

% time array
day = 3600*24; % seconds
year = 365*day; % seconds
dt = day; % time step day
tmax = year*1000*5; % max time ka
t = 0:dt:tmax; % time array

imax = length(t);

% initial bedrock
zb =ones(size(x)); % flat bedrock
E = .1/365; %m/year

% Setup initial temperature profile
T = ones(size(z))*Tsnaught + ((Qm*z)/k); %initial T is known based on geotherm
T(1) = T(1) + Tstep; % Surface Temperature after change, C

% plot animation
n=100; %number of plots
tplot = tmax/n; % time step of plot

%% run
for i = 2:imax % Calculates each time step
    
    dTdz = diff(T)/dz;
    Q = -k*dTdz;
    dQdz = diff(Q)/dz;
    dQdz = [0 dQdz 0];
    T = T - (1/(rho*c))*dQdz*dt;
    zb = zb+E;

% Plot solution
if (rem(t(i),tplot)==0) % decrease plotting frequency - speed up animation
figure(1), clf
plot(x,zb);
xlabel('Distance [m]')
ylabel('depth [m]')
set(gca,'YDIR','reverse')
set(gca,'fontsize',14,'fontname','arial')
title(['Erosion after ',num2str(t(i)/year),' years'])
set(gca,'YDIR','reverse')
axis([0 xmax 0 1000]) % hold axes constant
pause(0.05)
end

if (rem(t(i),tplot)==0) % decrease plotting frequency - speed up animation
figure(2), clf
plot(T,z);
hold on
plot(T,zb)
xlabel('Temperature [^oC]')
ylabel('depth [m]')
set(gca,'YDIR','reverse')
set(gca,'fontsize',14,'fontname','arial')
title(['Temperature evolution after ',num2str(t(i)/year),' years'])
set(gca,'YDIR','reverse')
axis([0 60 0 zmax]) % hold axes constant
pause(0.05)
end

end
