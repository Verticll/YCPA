winstyle = 'docked';
% winstyle = 'normal';

set(0,'DefaultFigureWindowStyle',winstyle)
set(0,'defaultaxesfontsize',18)
set(0,'defaultaxesfontname','Times New Roman')
% set(0,'defaultfigurecolor',[1 1 1])

% clear VARIABLES;
clear
global spatialFactor;
global c_eps_0 c_mu_0 c_c c_eta_0
global simulationStopTimes;
global AsymForcing
global dels
global SurfHxLeft SurfHyLeft SurfEzLeft SurfHxRight SurfHyRight SurfEzRight

% This is a 2D simulation of a traveling EM wave using FDTD

dels = 0.75;
spatialFactor = 1;

c_c = 299792458;                  % speed of light
c_eps_0 = 8.8542149e-12;          % vacuum permittivity
c_mu_0 = 1.2566370614e-6;         % vacuum permeability
c_eta_0 = sqrt(c_mu_0/c_eps_0);


tSim = 200e-15 % total time of the simulation
f = 130e12; % frequency
lambda = c_c/f;

xMax{1} = 20e-6;
nx{1} = 200; %Setting Size of space in x
ny{1} = 0.75*nx{1};  %Setting size of space in y to be 75% of X


Reg.n = 1;

mu{1} = ones(nx{1},ny{1})*c_mu_0;

epi{1} = ones(nx{1},ny{1})*c_eps_0;
%The line below is the inclusion and when removed removes the box and
%allows to wave to travel without interference

%Added assymetric pattern which I thought looked neat
epi{1}(125:150,55:95)= c_eps_0*11.3;
epi{1}(85:110,90:110)= c_eps_0*11.3;
epi{1}(85:150,60:80)= c_eps_0*11.3;



sigma{1} = zeros(nx{1},ny{1});
sigmaH{1} = zeros(nx{1},ny{1});

dx = xMax{1}/nx{1}; %get the space step
dt = 0.25*dx/c_c;   %get the time step
nSteps = round(tSim/dt*2); % get the amount of time steps in the simulation
yMax = ny{1}*dx; 
nsteps_lamda = lambda/dx %get the amount of wavelength steps in the simulation

movie = 1;
Plot.off = 0;
Plot.pl = 0;
Plot.ori = '13';
Plot.N = 100;
Plot.MaxEz = 1.1;
Plot.MaxH = Plot.MaxEz/c_eta_0;
Plot.pv = [0 0 90];
Plot.reglim = [0 xMax{1} 0 yMax];

%Setting up boundary conditions and source info
bc{1}.NumS = 1;
bc{1}.s(1).xpos = nx{1}/(4) + 1;
bc{1}.s(1).type = 'ss';
bc{1}.s(1).fct = @PlaneWaveBC;
% mag = -1/c_eta_0;
mag = 2; % Magnitude
phi = 0;    %phase
omega = f*2*pi;
betap = 0; %unsure
t0 = 30e-15;
%made the grating clearer
st = -0.05;
s = 0; % Reflection?
y0 = yMax/2;
sty = 1.5*lambda;
bc{1}.s(1).paras = {mag,phi,omega,betap,t0,st,s,y0,sty,'s'};

Plot.y0 = round(y0/dx);

bc{1}.xm.type = 'a';
bc{1}.xp.type = 'a';
bc{1}.ym.type = 'a';
bc{1}.yp.type = 'a';

pml.width = 20 * spatialFactor;
pml.m = 3.5;

Reg.n  = 1;
Reg.xoff{1} = 0;
Reg.yoff{1} = 0;

RunYeeReg % Runs the simulation and plotting






