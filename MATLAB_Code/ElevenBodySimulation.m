clear
clc
close all

%% Constants (All constants are in SI Units)

N = 11; % Number of bodies
G = 6.6743e-11; % Gravitational constant
m = [1988500e24,6.4171e23,48.685e23,5.97219e24,7.349e22,6.39e23,1.89813e27,5.683e26,8.681e25,1.024e26,1.30900e22]*1000'; % masses of the objects
n_dts_per_ode = 20; % Number of timesteps to solve with the ode solver during one iteration of the loop
R = 20; % Number of times to run RK4
dt = 6; % Timesteps
n_dts = 30; % Number of timesteps for one period of orbit
maxtime = dt*n_dts_per_ode; % Time span for one ode call
time_span=0:dt:maxtime; %Span ode will operate on
max_iter_to_plot=100; % How many points to trim the trajectory plot to
%% Initial conditions    
r=[0 0 0 0 0 0
  1.835679060203390E+07 -6.477899324095604E+07 -6.977492807261240E+06 3.710210432823309E+01 1.576057209541615E+01 -2.115250957994619E+00
  -1.887637152638366E+05 -1.087453135732388E+08 -1.482039133512773E+06 3.478559988580896E+01 -1.895996072366296E-01 -2.009817292668111E+00
  5.663095182784107E+07 1.362590372611622E+08 -6.947178111299872E+03 -2.799779113982546E+01 1.130984210949225E+01 -7.314488611189773E-04
  5.695003603541929E+07 1.360712632871835E+08 -3.945215585505217E+04 -2.742509968994457E+01 1.219383426385233E+01 2.293954446209234E-02
  7.226140454726927E+07 2.161676301702433E+08 2.757949164332956E+06 -2.206118545854885E+01 9.743131019701204E+00 7.453489395585069E-01
  7.309385491412139E+08 1.191406485942078E+08 -1.684833812295864E+07 -2.255559923496945E+00 1.352370922660542E+01 -5.638108121003071E-03
  1.205173185532243E+09 -8.455353900080341E+08 -3.325967557790118E+07 5.011387204192125E+00 7.901420062279396E+00 -3.366930787829872E-01
  2.014001339620603E+09 2.146255859985742E+09 -1.813384333001649E+07 -5.018439642050031E+00 4.356612663152282E+00 8.085950310893608E-02
  4.450760552111814E+09 -4.553229099109694E+08 -9.319709403937596E+07 5.161770702345839E-01 5.454981447847490E+00 -1.238582115111053E-01
  2.406342567477747E+09 -4.587627672010628E+09 -2.047987974797966E+08 4.946008581210541E+00 1.341889246804544E+00 -1.594378960221578E+00];


%% Plot Setup
myplot=figure('Position',[100 100 850 850]);hold on;
xlabel('X');ylabel('Y');zlabel('Z');
title('Model of Solar System')
daspect([1 1 1]); 
view(3);

% Plot of each orbit
for i=1:N
   body_traj(i)=plot3(r(i,1), r(i,2), r(i,3),'-'); hold on       % Orbit line
   body(i)=plot3(r(i,1), r(i,2), r(i,3),'o', 'MarkerSize',8);    % Planet Marker
end

t=0; % Absolute time
% whitebg('black');

grid on;
drawnow   
options = odeset('reltol',1e-3,'abstol',5e-2); % tolerance for ode solver


row_count = 0;
History = zeros(length(time_span),66);
XX = zeros(length(time_span),11);
YY = zeros(length(time_span),11);
ZZ = zeros(length(time_span),11);
VX = zeros(length(time_span),11);
VY = zeros(length(time_span),11);
VZ = zeros(length(time_span),11);
for i = 1:R
    [TTT,XXX] = ode89(@ElevenBody,time_span,r,options);
    steps=length(TTT);
    History=XXX(:,:);
    Size = size(XXX,1);

% Get new coordinates from solution    
    i=1:N;
    XX(:,i)=History(:,i);
    YY(:,i)=History(:,i+N);
    ZZ(:,i)=History(:,i+2*N);
    
    VX(:,i)= History(:,i+3*N);
    VY(:,i)= History(:,i+4*N);
    VZ(:,i)= History(:,i+5*N);
    
    t=t+dt;

% Replot points
    for ii=1:N
            body_traj(ii).XData=[body_traj(ii).XData XX(:,ii)'];
            body_traj(ii).YData=[body_traj(ii).YData YY(:,ii)']; %    
            body_traj(ii).ZData=[body_traj(ii).ZData ZZ(:,ii)']; % Update plots    
            set(body(ii), 'XData',XX(end,ii),'YData',YY(end,ii),'ZData',ZZ(end,ii));
    end

r=reshape(XXX(end,:),[11,6]); % Use current state as initial conditions for next portion of orbit

t=t+maxtime-dt; 
drawnow
end