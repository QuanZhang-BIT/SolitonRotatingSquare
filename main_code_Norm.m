% ---------- main function  -----------
%warning('off')
clear all
close all
clc
global a theta0 h Hh r Vol vol
global Nx Ny N damp Ndamp
global M J alfa
global k_theta k_s k_l 
global K_theta K_s
global T_out theta_out Ux_out Uy_out
global theta_st ex_st ey_st
global mu0 m Mm B Bm
global fai cc
global phase_shift

% This code simulate the response of Nx by Ny squares connected
%% User defined parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% structural parameters
a = 10*10^(-3); % lattice constant
theta0 = 5/180*pi; % initial rotate angle
h = 0.002; % out-of-plane thickness of the lattice
Hh = h/a; % normalized out-of-plane thickness
r = a/3/cos(theta0)/sqrt(2); % radius of the magnetic core
Vol = pi*r^2*h; % volume of the magnetic core
vol = Vol/a^3; % normalized volume
Nx = 60; % number of units in x direction
Ny = 120; % number of units in y direction
N = Nx*Ny; % number of units in the lattice
damp = 0.2; % damping factor
Ndamp = 15; % number of units with damping near the boundary
%Nexci = 10; % number of units being excited 
M = 74.3*10^(-3)*h; % mass of the square
J = 599.8*10^(-9)*h; % moment of inertia of the square
alfa = a*sqrt(M/J); % normalized moment of inertia
k_l = 186.4*10^3*h; % longitudinal stiffness of the hinge
k_s = 48.0*10^3*h; % shear stiffness of the hinge
K_s = k_s/k_l; % normalized shear stiffness
k_theta = 12.02*10^(-3)*h; % rotate stiffness of the hinge
K_theta = k_theta/(k_l*a^2); % normalized rotate stiffness

% magnetic parameter
mu0 = 4*pi*10^(-7); % magnetic permeability
m = 64*10^3; % magnetization of the magnetic core
Mm = m*sqrt(mu0*a/k_l); % normalized magnetization
B = 0*10^(-3); % strength of the appled magnetic field
Bm = B/sqrt(mu0*k_l/a); % normalized magnetic field

% wave parameters
fai = 0/180*pi; % angle for the propagation direction
cc = 11.8; % wave velocity

% simulation time & animation play rate
phase_shift = 25*a;
Time = (phase_shift+(Nx-10)*a)/cc*sqrt(k_l/M); %duration of the simulatios
Framenumber = Time/20; %number of frames used to generate the movies

%% magnetic field_induced predeformation
[theta_st,ex_st,ey_st] = magnetic_predeformation(B);

%% obtain theoretical solution
[T_out,theta_out,Ux_out,Uy_out] = input_generator_Norm(cc,B);

%% initial condition
% y0 = 0.0000001.*rand(4*N,1);
y0 = 0.*rand(6*N,1);

%% integrate the system of ODEs with ODE45
tic
options=odeset('RelTol',1e-5,'AbsTol',1e-7,'Stats','off');
[T,Y] = ode45('ODEsNorm',[0,Time], y0);
[Timelength, Distance]=size(Y);
timeelapse=toc;
disp(' '); disp(['Approximate running time is: ', num2str(timeelapse)])

%% post processing data
clear Ux Uy theta
X = 1:Nx;

for j=X
    Ux(:,j) = Y(:,6*((round(Ny/2)-1)*Nx+j-1)+1); % displacement along the center line of the lattice
    Uy(:,j) = Y(:,6*((round(Ny/2)-1)*Nx+j-1)+3);
    theta(:,j) = Y(:,6*((round(Ny/2)-1)*Nx+j-1)+5);
end
Rmin = min(theta_out);Rmax = max(theta_out);
Uxmin = min(Ux_out);Uxmax = max(Ux_out);
Uymin = min(Uy_out);Uymax = max(Uy_out);
%% save data
save theta;
save Ux;
save Uy;
save T;
save Y;

%% plot displacement and rotation signal gif
fig = figure(1);
clf
set(fig,'position', [0, 0, 1500,300])
set(gcf,'Color',[1,1,1])
top_folder = mkdir('results');
sol = strcat('B=',num2str(B),'_c_',num2str(cc),'_numNx_',num2str(Nx),'_numNy_',num2str(Ny),'_Norm','_BoundaryDamping');
mkdir('results/',sol);
file = [strcat('Ux_','B=',num2str(B)),'_fai_',num2str(fai),'_c_',num2str(cc),'_numNx_',num2str(Nx),'_numNy_',num2str(Ny),'_numNdamp_',num2str(Ndamp),'_damp_',num2str(damp)];
filename1 = ['results/',sol,'/',file,'.gif'];
count = 0;
for t=1:Timelength
    if T(t)*cc/sqrt(k_l/M)/a > count*1-1
        count = count + 1;
        hold off
        skip=1;
        plot([2:skip:Nx],Ux(t,2:skip:end),'-ko','markerfacecolor','r','MarkerSize',7);% descrete model
        hold on;
        plot((-T_out*cc/sqrt(k_l/M)+cc*T(t)/sqrt(k_l/M)+a)/a, Ux_out, '-r','linewidth',1.5); % theoretical solution
        ylabel('Disp., \itu_{\itx}/\ita','FontName','Arial','fontsize',28)
        ylim([Uxmin-0.1,Uxmax+0.05]);
%         ylabel('Angle, \it\theta','FontName','Arial','fontsize',28)
%         ylim([Rmin-0.04,+0.06]);
        xlabel('Position, \itx/\ita','FontName','Arial','fontsize',28)
        set(gca,'FontName','Arial','fontsize',24)
        xlim([-Nx,Nx])
        box on
        pause(0.001);
        frame = getframe(fig);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,1024);
        if count == 1
            imwrite((imind),cm,filename1,'gif','DelayTime',0.01, 'Loopcount',inf);
        else
            imwrite(uint8(imind),cm,filename1,'gif','DelayTime',0.01,'WriteMode','append');
        end
    end
end

%% plot displacement and rotation signal fig
fig = figure(1);
clf
set(fig,'position', [0, 0, 500,300])
set(gcf,'Color',[1,1,1])
T_target = [15 30 45]+phase_shift/a-1;
nT_target = length(T_target);

for i = 1:nT_target
    min_dT = min(abs(cc*T/sqrt(k_l/M)/a-T_target(i)));
    t = find(abs(cc*T/sqrt(k_l/M)/a-T_target(i))==min_dT);
    skip=1;
    hold on;
    plot((-T_out*cc/sqrt(k_l/M)+cc*T(t)/sqrt(k_l/M)+a)/a, theta_out, '-k','linewidth',1.5); % theoretical solution

    set(gca,'FontName','Arial','fontsize',24)
    xlim([0,Nx])
    box on
end
for i = 1:nT_target
    min_dT = min(abs(cc*T/sqrt(k_l/M)/a-T_target(i)));
    t = find(abs(cc*T/sqrt(k_l/M)/a-T_target(i))==min_dT);
    skip=1;
%     scatter([2:skip:Nx],Ux(t,2:skip:end),30,'markeredgecolor','k','markerfacecolor',[47/255 71/255 156/255],'MarkerFaceAlpha',0.7);% descrete model
%     ylim([Uxmin-0.1,0.4]);
%     %ylim([Uxmin-0.03,0.1]);
    
    scatter([2:skip:Nx],theta(t,2:skip:end),30,'markeredgecolor','k','markerfacecolor',[255/255 128/255 0/255],'MarkerFaceAlpha',0.6);% descrete model
    ylim([-0.05,0.31]);
    set(gca,'FontName','Arial','fontsize',24)
    xlim([0,Nx])
    box on
end



%% plot displacement field gif
load T;
load Y;
[Timelength, Distance]=size(Y);
disp_theta = Y(:,5:6:Distance);

x = [0 0 0 0];
y = [0 0 0 0];
c = 0;
fig = figure(2);
clf
set(fig,'position', [100, 100, 700,700])
set(gcf,'Color',[1,1,1])
axis equal
axis off

top_folder = mkdir('results');
sol = strcat('B=',num2str(B),'_c_',num2str(cc),'_numNx_',num2str(Nx),'_numNy_',num2str(Ny),'_Norm','_BoundaryDamping');
mkdir('results/',sol);
file = [strcat('UxField_','B=',num2str(B)),'_fai_',num2str(fai),'_c_',num2str(cc),'_numNx_',num2str(Nx),'_numNy_',num2str(Ny),'_numNdamp_',num2str(Ndamp),'_damp_',num2str(damp)];
filename1 = ['results/',sol,'/',file,'.gif'];
count1 = 0;
for t=1:Timelength
    if T(t)*cc/sqrt(k_l/M)/a > count1*1-1
        count1 = count1 + 1;
        hold off
        for i=1:Ny
            for j=1:Nx
                x(1) = a/2/cos(theta0)*cos(theta0)+(j-1)*a;
                y(1) = -a/2/cos(theta0)*(-1)^(i+j)*sin(theta0)+(i-1)*a;
                x(2) = a/2/cos(theta0)*(-1)^(i+j)*sin(theta0)+(j-1)*a;
                y(2) = a/2/cos(theta0)*cos(theta0)+(i-1)*a;
                x(3) = -a/2/cos(theta0)*cos(theta0)+(j-1)*a;
                y(3) = a/2/cos(theta0)*(-1)^(i+j)*sin(theta0)+(i-1)*a;
                x(4) = -a/2/cos(theta0)*(-1)^(i+j)*sin(theta0)+(j-1)*a;
                y(4) = -a/2/cos(theta0)*cos(theta0)+(i-1)*a;
                c = disp_theta(t,(i-1)*Nx+j)/max(max(theta));
                patch(x,y,c);
                caxis([0,1]);
                hold on;
            end
        end
        pause(0.001);
        frame = getframe(fig);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,1024);
        if count1 == 1
            imwrite((imind),cm,filename1,'gif','DelayTime',0.01, 'Loopcount',inf);
        else
            imwrite(uint8(imind),cm,filename1,'gif','DelayTime',0.01,'WriteMode','append');
        end
    end
end

%% plot displacement field fig
load T;
load Y;
[Timelength, Distance]=size(Y);
disp_theta = Y(:,3:6:Distance);

x = [0 0 0 0];
y = [0 0 0 0];
c = 0;
fig = figure(2);
clf
set(fig,'position', [100, 100, 700,700])
set(gcf,'Color',[1,1,1])
axis equal
axis off

T_target = [15 30 45]+phase_shift/a-1;
nT_target = length(T_target);
for ii =2
    min_dT = min(abs(cc*T/sqrt(k_l/M)/a-T_target(ii)));
    t = find(abs(cc*T/sqrt(k_l/M)/a-T_target(ii))==min_dT);
    skip=1;
    for i=1:Ny
        for j=1:Nx
            x(1) = a/2/cos(theta0)*cos(theta0)+(j-1)*a;
            y(1) = -a/2/cos(theta0)*(-1)^(i+j)*sin(theta0)+(i-1)*a;
            x(2) = a/2/cos(theta0)*(-1)^(i+j)*sin(theta0)+(j-1)*a;
            y(2) = a/2/cos(theta0)*cos(theta0)+(i-1)*a;
            x(3) = -a/2/cos(theta0)*cos(theta0)+(j-1)*a;
            y(3) = a/2/cos(theta0)*(-1)^(i+j)*sin(theta0)+(i-1)*a;
            x(4) = -a/2/cos(theta0)*(-1)^(i+j)*sin(theta0)+(j-1)*a;
            y(4) = -a/2/cos(theta0)*cos(theta0)+(i-1)*a;
            v = [x(1) y(1);x(2) y(2);x(3) y(3);x(4) y(4)];
            f = [1 2 3 4];
            c = disp_theta(t,(i-1)*Nx+j)/max(max(Ux));
%             patch(x,y,c);
            patch('Faces',f,'Vertices',v,'FaceVertexCData',c,'FaceColor','flat','EdgeColor',[0 0 0]/256,'LineWidth',0.2);
            caxis([-1,1]);
            hold on;
        end
    end
end

%% plot spatial-temporal map

fig = figure(2);
clf
set(fig,'position', [100, 100, 700,300])
t = cc*T(1:5:end)/sqrt(k_l/M)/a-phase_shift/a+1;
x = 1:60;
[tX,tT] = meshgrid(x,t);
pcolor(tX,tT,Uy(1:5:end,:));
caxis([-0.27,0.27]);
shading interp;
xlim([2,60]);
ylim([0,50]);
set(gca,'XTick',20:20:60);
set(gca,'XTicklabel',{'20','40','60'});
set(gca,'YTick',0:15:50);
set(gca,'YTicklabel',{'0','15','30','45'});
set(gca,'FontName','Arial','fontsize',32)


% close all
