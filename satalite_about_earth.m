%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg



% this program is used for showing a satalite about the earth
close all; clear all; clc;
%% constants
G=6.67428E-11; % gravitational constant
hz=15000; % simulation frequancy
sim_hz=20;  % simulation step
%% Intial condition
M=[5.972*10^24,1000];       % [M1,M2]
R1=[0;0;0];                 % position of M(1)
R2=[8000e3;0;0];                 % position of M(2)
I=[0,0,0];              % location of initial axis
V1=[0;0;0];                 % velocity of M(1)
V2=[0;6;3]*1e3;                 % velocity of M(2)
% image (inter your URL of the image)
image_file = 'D:\4th year of Aerospace\1st\Orbital Mechanics\AER-427, Orbital Mechanics, Mohamed Mohamed Elsayed,SCE 2, BN 13  By MATLAB\week 11\satalite about earth with tracking/03.jpg';
%% RK4 parameter
tf=3650*5;   % final time of soution
dt=2;            % time step
X0=[R1;R2;V1;V2];
B=[0;0;0;0;0;0;0;0;0;0;0;0];
sol(1:12,1)=X0;
order=12;
%% solution by RK4
for n=1:length(0:dt:tf)
    b=G*M(2)/(norm(sol(1:3,n)-sol(4:6,n)))^3;
    c=-G*M(1)/(norm(sol(1:3,n)-sol(4:6,n)))^3;
    A=[0,0,0,0,0,0,1,0,0,0,0,0; ...
        0,0,0,0,0,0,0,1,0,0,0,0; ...
        0,0,0,0,0,0,0,0,1,0,0,0; ...
        0,0,0,0,0,0,0,0,0,1,0,0; ...
        0,0,0,0,0,0,0,0,0,0,1,0; ...
        0,0,0,0,0,0,0,0,0,0,0,1;...
        -b,0,0,b,0,0,0,0,0,0,0,0; ...
        0,-b,0,0,b,0,0,0,0,0,0,0; ...
        0,0,-b,0,0,b,0,0,0,0,0,0; ...
        -c,0,0,c,0,0,0,0,0,0,0,0; ...
        0,-c,0,0,c,0,0,0,0,0,0,0; ...
        0,0,-c,0,0,c,0,0,0,0,0,0 ];
    [ XX ] = RK4( A,B,sol(1:12,n),dt,n*dt,(n+1)*dt,order );
    sol(1:12,n+1)=XX(1:12,2);
end
R1_x=sol(1,:);
R1_y=sol(2,:);
R1_z=sol(3,:);
R2_x=sol(4,:);
R2_y=sol(5,:);
R2_z=sol(6,:);
V1_x=sol(7,:);
V1_y=sol(8,:);
V1_z=sol(9,:);
V2_x=sol(10,:);
V2_y=sol(11,:);
V2_z=sol(12,:);
%% center of masses parameters
r=[R2_x;R2_y;R2_z]-[R1_x;R1_y;R1_z];                                            % the distance betweem M1 & M2
Rc=(M(1)*[R1_x;R1_y;R1_z]+M(2)*[R2_x;R2_y;R2_z])/sum(M);                        % location of center of masses
Vc=(M(1)*[V1_x;V1_y;V1_z]+M(2)*[V2_x;V2_y;V2_z])/sum(M);                        % Vc = constant
Ac=(M(1)*G*M(2)*r.^3.*r-M(1)*G*M(2)*r.^3.*r)/sum(M);    % for check Ac = 0
for H=1:length(Rc(1,:))
    acc1(1:3,H)=(-1)^(1)*G*M(1)/(norm([R2_x(H);R2_y(H);R2_z(H)]'-[R1_x(H);R1_y(H);R1_z(H)]'))^3*([R2_x(H);R2_y(H);R2_z(H)]-[R1_x(H);R1_y(H);R1_z(H)]); % acceleration of M1
    acc2(1:3,H)=(-1)^(2)*G*M(2)/(norm([R2_x(H);R2_y(H);R2_z(H)]'-[R1_x(H);R1_y(H);R1_z(H)]'))^3*([R2_x(H);R2_y(H);R2_z(H)]-[R1_x(H);R1_y(H);R1_z(H)]); % acceleration of M2
end
%% velocity and acceleration magnituides
for h=1:length(Rc(1,:))
        MagV1(1,h)=norm([V1_x;V1_y;V1_z]);  % velocity magnituide of M1
        MagV2(1,h)=norm([V2_x;V2_y;V2_z]);  % velocity magnituide of M2
        MagA1(1,h)=norm(acc1(1:3,h));  % acceleration magnituide of M1
        MagA2(1,h)=norm(acc2(1:3,h));  % acceleration magnituide of M2
end
%% projected path of the satalite on the earth
erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)
for P=1:length(Rc(1,:))
    points = intersectLineSphere([0,0,0,R2_x(P)-R1_x(P)+(P-1)*erot*dt,R2_y(P)-R1_y(P)+(P-1)*erot*dt,R2_z(P)-R1_z(P)+(P-1)*erot*dt], [0,0,0,6400e3]);
    xx(P)=points(2,1);
    yy(P)=points(2,2);
    zz(P)=points(2,3);
end
%% plotting
%--------------------------------------------------------------------------------------------------------------------------------------------------------
figure(1);
for p=1:sim_hz:length(R1_x)
% for p=1
    % Options
    space_color = 'k';
    npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
    alpha   = 1; % globe transparency level, 1 = opaque, through 0 = invisible
    % Earth texture image
    % Anything imread() will handle, but needs to be a 2:1 unprojected globe
    % Mean spherical earth
    erad    = 6371008.7714; % equatorial radius (meters)
    prad    = 6371008.7714; % polar radius (meters)
    %GMST0 = []; % Don't set up rotatable globe (ECEF)
    GMST0 = 4.89496121282306 + (p-1)*erot*dt; % Set up a rotatable globe at J2000.0
    set(gcf,'Color','w');
    % M2 trajectory
    plot3(R2_x-R1_x,R2_y-R1_y,R2_z-R1_z,'color','g','LineWidth',2);
    hold on;
    grid on;
    xlabel('X','Fontsize',18);
    ylabel('Y','Fontsize',18);
    zlabel('Z','Fontsize',18);
    title('Satalite about Earth','Fontsize',18);
    xlim auto;
    ylim auto;
    zlim auto;
    view(-84,44);
    % M2 start
    plot3(R2_x(p)-R1_x(p),R2_y(p)-R1_y(p),R2_z(p)-R1_z(p),'o','color','cyan','LineWidth',5);
    % Create wireframe globe
    % Create a 3D meshgrid of the sphere points using the ellipsoid function
    [x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
    globe = surf(x+R1_x(p), y+R1_y(p), -z+R1_z(p), 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
    if ~isempty(GMST0)
        hgx = hgtransform;
        set(hgx,'Matrix', makehgtform('zrotate',GMST0));
        set(globe,'Parent',hgx);
    end
    % Texturemap the globe
    % Load Earth image for texture map
    cdata = imread(image_file);
    % Set image as color data (cdata) property, and set face color to indicate
    % a texturemap, which Matlab expects to be in cdata. Turn off the mesh edges.
    set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
    %Projected path of the satalite on the earth
    plot3(xx,yy,zz,'r','LineWidth',2);
    legend('Satalite trajectory','Satalite','Projected path of the satalite on the earth');
    pause(1/hz);
    hold off;
end
%--------------------------------------------------------------------------------------------------------------------------------------------------------