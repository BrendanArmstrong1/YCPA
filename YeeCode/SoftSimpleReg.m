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


%% B
% This file is simulating a wave going through transparaent medium, there
% is an object in the medium with a different refractive index that the
% wave interacts with.





dels = 0.75;
spatialFactor = 1;

c_c = 299792458;                  % speed of light
c_eps_0 = 8.8542149e-12;          % vacuum permittivity
c_mu_0 = 1.2566370614e-6;         % vacuum permeability
c_eta_0 = sqrt(c_mu_0/c_eps_0);


tSim = 200e-15   %sumilation duration
f = 230e12;      % frequency of the wave
lambda = c_c/f;  % wavelength

xMax{1} = 20e-6;  % distance the wave travels
nx{1} = 200;      % number of divisions in X axis
ny{1} = 0.75*nx{1}; % number of divisions in Y axis


Reg.n = 1;

mu{1} = ones(nx{1},ny{1})*c_mu_0;  % permeability of material throughout the grid

epi{1} = ones(nx{1},ny{1})*c_eps_0; % permittivity throughout the grid
epi{1}(125:150,55:95)= c_eps_0*11.3; % This adds the inclusions. commenting
                                     % out does work
                                     
%adding more inclusions
epi{1}(100:101,55:95)= c_eps_0*11.3;
epi{1}(95:96,55:95)= c_eps_0*11.3;
epi{1}(92:93,55:95)= c_eps_0*11.3;
                                     
                                     
                                     
sigma{1} = zeros(nx{1},ny{1});
sigmaH{1} = zeros(nx{1},ny{1});

dx = xMax{1}/nx{1};  % division size
dt = 0.25*dx/c_c;    % time size
nSteps = round(tSim/dt*2);
yMax = ny{1}*dx;
nsteps_lamda = lambda/dx

movie = 1;
Plot.off = 0;
Plot.pl = 0;
Plot.ori = '13';
Plot.N = 100;
Plot.MaxEz = 1.1;
Plot.MaxH = Plot.MaxEz/c_eta_0;
Plot.pv = [0 0 90];
Plot.reglim = [0 xMax{1} 0 yMax];

%% C
% ii) The bc structure seems to be setting boundary conditions for the yee
% cells.


bc{1}.NumS = 2; % 2 sources
bc{1}.s(1).xpos = nx{1}/(4) + 1; % sets up the start of the wave visibility
                                 % in this case, the incoming wave is only
                                 % visible 1 quarter the x distance.
                                 
bc{1}.s(1).type = 'ss';
bc{1}.s(1).fct = @PlaneWaveBC;

bc{1}.s(2).xpos = nx{1}/(4) + 1;                                 
bc{1}.s(2).type = 'ss';
bc{1}.s(2).fct = @PlaneWaveBC;
% mag = -1/c_eta_0;
mag = 1; % magnitude
phi = 0; % seems to effect permeability
omega = f*2*pi; %wavelength
betap = 20;
t0 = 30e-15; % starting distance of the pulse. larger = further away
st = -.05; % duration of the pulse. setting to -.05 creates the wave ahead of time
s = 0;
y0 = yMax/2;  % altitude of the wave

sty = 1.5*lambda; % width of the wave
bc{1}.s(1).paras = {mag,phi,omega,betap,t0,st,s,y0,sty,'s'}; % implements the settings.
bc{1}.s(2).paras = {mag,phi,omega,betap,2*t0,st,s,y0,sty,'s'}; % second source delayed from first one

Plot.y0 = round(y0/dx);

% These settings seem to adjust how the boundaries work, weather they allow
% waves to pass through, wrap around, or reflect.
bc{1}.xm.type = 'a';
bc{1}.xp.type = 'a';
bc{1}.ym.type = 'a';
bc{1}.yp.type = 'a';



pml.width = 20 * spatialFactor;
pml.m = 3.5;

Reg.n  = 1;
Reg.xoff{1} = 0;
Reg.yoff{1} = 0;

RunYeeReg






