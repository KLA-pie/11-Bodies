clear all;clc;close all;
global G m
record=0; % shall we record video?



%% constants

N=11;
G=6.6743e-11; % gravitational constant
m=[1988500e24,6.4171e23,48.685e23,5.97219e24,7.349e22,6.39e23,1.89813e27,5.683e26,8.681e25,1.024e26,1.30900e22]'; % masses of the objects
n_dts_per_ode = 2; % number of timesteps to solve with ode solver during one iteration of while time loop
max_iter_to_plot=100; % how many points to trim the trajectory plot to
%% Initial conditions
% arranging objects -- optional, you can arrange them as you like
R=3; % radius of the circumference object will sit on initially

i=1:N; 
thetas=(360*i/N)';   % polar angles fo each object's initial position

X = R*cosd(thetas);         % Cartesian coordinates of the objects
Y = R*sind(thetas) ;        
Z = 0*thetas;               
        
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


% Finding (configuration) Potential energy of each mass in the gravitational field due to all other masses
PE=zeros(N,1);
for i=1:N
    for j=1:N
        if i~=j
            PE(i)=PE(i)+G*m(i)./norm(r(j,:)-r(i,:)); % P.E. of mass i due to all j's
        end
    end
end
v0=sqrt(PE/2); % orbital Keplerian speeds for circular orbit (PE=2*KE) for each mass

% velocity components
vx = -v0.*sind(thetas);
vy = v0.*cosd(thetas);
vz = v0*0; % start motion in one plane

inits=r; % initial conditions for each mass -- all in 1 column, 6 rows for ode solvers

%%%%%%% end initial conditions %%%%%%%%%%%


%% COM (barycentre)

COM=sum(r(1:N,1:3).*m(1:N),1)/sum(m); % location of barycentre

%% setup plotting
myplot=figure('Position',[100 100 850 850]);hold on; % figure window
COM_plot=plot3(COM(:,1),COM(:,2),COM(:,3),'*'); hold on;
% axis auto
xlabel('X');ylabel('Y');zlabel('Z');
daspect([1 1 1]); 
view(3); % change to view(3) for 3D projection
n_dts=360*1; % how many of timesteps per one period of the orbit

dt=1 % timestep: n_dts steps per period of period

% plots for each orbit of each mass
for i=1:N
   body_traj(i)=plot3(X(i), Y(i),Z(i),'-'); hold on       % orbit line
   body(i)=plot3(X(i), Y(i), Z(i),'o', 'MarkerSize',8);    % mass marker
end

t=0; % abs time
whitebg('black');

grid on;
drawnow   

maxtime=dt*n_dts_per_ode; % time span for one ode45 call

time_span=0:dt:maxtime;
options = odeset('reltol',1e-6,'abstol',1e-4); % tolerance for ode solver


n=0; % iteration count


while isgraphics(myplot) % infinite DKD leapfrog time loop (while plot open)
    x_history=[];
    [TTT,XXX] = ode(@ElevenBody,time_span,inits,options);
    steps=length(TTT);
    x_history=XXX(:,:);

% get new coordinates from solution    
    i=1:N;
    XX(:,i)=x_history(:,i);
    YY(:,i)=x_history(:,i+N);
    ZZ(:,i)=x_history(:,i+2*N);
    
    VX(:,i)= x_history(:,i+3*N);
    VY(:,i)= x_history(:,i+4*N);
    VZ(:,i)= x_history(:,i+5*N);
    
    n=n+1; % number of total ode runs performed
    t=t+dt;

% redrawind all plots
    for ii=1:N
            body_traj(ii).XData=[body_traj(ii).XData XX(:,ii)'];
            body_traj(ii).YData=[body_traj(ii).YData YY(:,ii)']; %    
            body_traj(ii).ZData=[body_traj(ii).ZData ZZ(:,ii)']; % update plots    
            set(body(ii), 'XData',XX(end,ii),'YData',YY(end,ii),'ZData',ZZ(end,ii)); % update mass' markers
            
        if n*steps>=max_iter_to_plot % trimming the trajectory line to max_iter_to_plot
            body_traj(ii).XData=body_traj(ii).XData(end-max_iter_to_plot+1:end);
            body_traj(ii).YData=body_traj(ii).YData(end-max_iter_to_plot+1:end);
            body_traj(ii).ZData=body_traj(ii).ZData(end-max_iter_to_plot+1:end);
        end
    end  %%%%%%%%%%%%%%%%% loop ii   ends
  
r=[XX(end,:)' YY(end,:)' ZZ(end,:)'];
COM=sum(r(i,1:3).*m(i),1)/sum(m); % new location of barycentre
set(COM_plot,'XData', COM(end,1),'YData', COM(end,2),'ZData', COM(end,3)); 

inits=XXX(end,:)'; % use current state as initial conditions for next portion of orbit

t=t+maxtime-dt; 
time.String=['t = ' num2str(round(t,3))];
drawnow

if record==1
    A=getframe(myplot);
    writeVideo(vid,A);
end


%% energy
% Kinetic Energy
V0=sqrt(VX(end,:).^2+VY(end,:).^2+VZ(end,:).^2);
KE=sum(1/2*m'.*V0.^2,2); % total Kintetic Energy


% Potential Energy
PE=zeros(N,1);
for i=1:N
    for j=1:N
        if i~=j
            PE(i)=PE(i)-G*m(i)./norm(r(j,:)-r(i,:)); % P.E. of mass i due to all j's
        end
    end
end
PE_tot=sum(PE);
E_txt.String=['Total Energy = ' num2str(KE+PE_tot)];

 end % while time loop ends

if record==1 
    close(vid); % close video file
end 
