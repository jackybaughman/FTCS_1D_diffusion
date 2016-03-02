% Geothermal gradient - relaxing after major erosion event
% 
% FTCS
% 
% JSB code update 
clear all
figure(2)
clf

%% Constants
Ts = 20; %initial surface T C

k = 2; % conductivity
rho = 2000; % density, kg/m^3
c = 2000; % Specific heat capacity, J/kg*K
kappa = k/(rho*c); % Thermal diffusivity, m^2/s

Qm = .04; % mantle heat flux

%% Arrays
% z array
dz = .01; % z step m
zmax = 200; % zmax m
z = 0:dz:zmax;% z array

% time array
day = 3600*24; % seconds
year = 365*day; % seconds
dt = day; % time step day
tmax = year*1000*1; % max time ka
t = 0:dt:tmax; % time array

imax = length(t);

% Setup initial temperature profile
T = ones(size(z))*Ts + ((Qm*z)/k); %initial T is known based on geotherm
T(1) = Ts;

% initial bedrock
Edot = -.000001/365; %m/year

% plot animation
n=100; % number of plots
tplot = tmax/n; % time step of plot

%% run
for i = 1:imax % Calculates each time step
    
    dTdz = diff(T)/dz;
    Q = -k*dTdz;
    dQdz = diff(Q)/dz;
    dQdz = [0 dQdz 0];
    T = T - (1/(rho*c))*dQdz*dt;
    znew = z+(Edot*dt);
    T = interp1(z,T,znew);
    z = znew;
    T(1) = Ts;

% Plot solution
if (rem(t(i),tplot)==0) % decrease plotting frequency - speed up animation
figure(2), clf
plot(T,z);
hold off
xlabel('Temperature [^oC]')
ylabel('depth [m]')
set(gca,'fontsize',14,'fontname','arial')
title(['Temperature evolution after ',num2str(t(i)/year),' years'])
set(gca,'YDIR','reverse')
axis([0 60 -100 200]) % hold axes constant
pause(0.05)
end

end

