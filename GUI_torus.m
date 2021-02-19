clear;
close all; 
clc


%% Initialize window

set(0,'Units','pixels')
dim = get(0,'ScreenSize');
fig_1 = [1];

fig_1 = figure('doublebuffer','on','Position',[0,35,dim(3)-200,dim(4)-100],...
            'Name',' Torus and wall intersection',...
            'NumberTitle','off');
ha3d = axes('Parent', fig_1,'Units', 'normalized','Position', [0.1 0.1 .8 .8]);
%set(fig_1, 'CurrentAxes', ha3d)
setappdata(0,'ha3d',[ha3d])

hold on;
%light('Position',[-1 0 0]);
light                               % add a default light
% daspect([1 1 1])                    % Setting the aspect ratio
view(135,25)
xlabel('X'),ylabel('Y'),zlabel('Z');
title ('Torus and plane');
axis([-8 8 -8 8 -8 8]);

plot3([-8,8],[-8,-8],[-8,-8],'k')
grid
%%
R = 2.5; % outer radius of torus
r = 1; % inner tube radius

% the plane - in our case wall
alpha = 0;      % rotation around Z axis, azimuth angle
phi = 0;    % roatation around Y axis, elevation angle
rho = 2;    % distance of the plane from the center of the torus

setappdata(0,'old_alpha',[alpha])
setappdata(0,'old_phi',[phi])
setappdata(0,'old_rho',[rho])

u_vector = [rho *  cosd(alpha) * cosd(phi); rho* sind(alpha) * cosd(phi); rho* sind(phi)] % vector 

Q = [rho * cosd(alpha) * cosd(phi) rho*sind(alpha)*cosd(phi) rho*sind(phi)]
p = Q +[sind(alpha) -cosd(alpha) 0]
p2 = (Q +[-cosd(alpha)*sind(phi)  -sind(alpha)*sind(phi) cosd(phi)])

t_axis = (p - Q)'
w_axis = (p2 - Q)'

setappdata(0,'old_t_axis',[t_axis])
setappdata(0,'old_w_axis',[w_axis])

% equation of plane
% helper equations
x_q = rho * cosd(alpha) * cosd(phi);
y_q = rho * sind(alpha) * cosd(phi);
z_q = rho * sind(phi);
setappdata(0,'old_x_q',[x_q])
setappdata(0,'old_y_q',[y_q])
setappdata(0,'old_z_q',[z_q])

x = x_q + t_axis * sind(alpha) - w_axis * cosd(alpha) * sind(phi)
y = y_q - t_axis * cosd(alpha) - w_axis * sind(alpha) * sind(phi)
z = z_q + w_axis * cosd(phi)

%% design plane
setappdata(0,'old_x',[x])
setappdata(0,'old_y',[y])
setappdata(0,'old_z',[z])
design_plane(x, y, z)

%% torus
setappdata(0,'old_R',[R])
setappdata(0,'old_r',[r])

tor(R, r)
hold on

%% Parameter Panel
param = uipanel(fig_1','units','normalized', 'Position',[0.0156 0.05625 0.245 ,0.31], ...
    'Title','Parameters','FontSize',12);
setappdata(0,'param',[param])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_slider = uicontrol(param,'style','slider',...
    'Max',10,'Min',0,'Value',2.5,...
    'SliderStep',[0.05 0.2],...        
    'callback',@R_slider_button_press,...
    'Position',[105 156 120 18]);
R_min = uicontrol(param,'style','text',...
    'String','0',...
    'Position',[75 157 25 14]); 
R_max = uicontrol(param,'style','text',...
    'String','+5',...
    'Position',[230 157 25 14]); 
R_text = uicontrol(param,'Style','text',...  
    'String','R - Outer radius of torus',...                
    'Position',[5 156 70 25]); 
R_edit = uicontrol(param,'style','edit',...
    'String',R_slider.Value,...
    'Position',[265 154 30 25]); 
setappdata(0,'R_slider',[R_slider])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_slider = uicontrol(param,'style','slider',...
    'Max',10,'Min',0,'Value',1,...
    'SliderStep',[0.05 0.2],...         
    'callback',@r_slider_button_press,...
    'Position',[105 126 120 18]);
r_min = uicontrol(param,'style','text',...
    'String','0',...
    'Position',[75 127 25 14]); 
r_max = uicontrol(param,'style','text',...
    'String','+5',...
    'Position',[230 127 25 14]); 
r_text = uicontrol(param,'Style','text',...  
    'String','r - inner tube radius',...      
    'Position',[5 126 70 25]); 
r_edit = uicontrol(param,'style','edit',...
    'String',r_slider.Value,...
    'Position',[265 124 30 25]); 
setappdata(0,'r_slider',[r_slider])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_slider = uicontrol(param,'style','slider',...
    'Max',360,'Min',0,'Value',0,...
    'SliderStep',[0.00278 0.02],...        
    'callback',@alpha_slider_button_press,...
    'Position',[105 96 120 18]);
alpha_min = uicontrol(param,'style','text',...
    'String','0',...
    'Position',[75 97 25 14]); 
alpha_max = uicontrol(param,'style','text',...
    'String','+360',...
    'Position',[230 97 25 14]); 
alpha_text = uicontrol(param,'Style','text',... 
    'String','a - rotation around Z axis',...                 
    'Position',[5 96 70 25]); 
alpha_edit = uicontrol(param,'style','edit',...
    'String', alpha_slider.Value,...
    'Position',[265 94 30 25]);
setappdata(0,'alpha_slider',[alpha_slider])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi_slider = uicontrol(param,'style','slider',...
    'Max',360,'Min',0,'Value',0,...
    'SliderStep',[0.00278 0.02],...       
    'callback',@phi_slider_button_press,...
    'Position',[105 66 120 18]);
phi_min = uicontrol(param,'style','text',...
    'String','0',...
    'Position',[75 67 25 14]);
phi_max = uicontrol(param,'style','text',...
    'String','+360',...
    'Position',[230 67 25 14]); 
phi_text = uicontrol(param,'Style','text',...  
    'String','phi - rotation around Y axis',...                 
    'Position',[5 66 70 25]); 
phi_edit = uicontrol(param,'style','edit',...
    'String',phi_slider.Value,...
    'Position',[265 64 30 25]);
setappdata(0,'phi_slider',[phi_slider])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho_slider = uicontrol(param,'style','slider',...
    'Max',10,'Min',-10,'Value',2,...
    'SliderStep',[0.005 0.2],...         
    'callback',@rho_slider_button_press,...
    'Position',[105 36 120 18]);
rho_min = uicontrol(param,'style','text',...
    'String','-10',...
    'Position',[75 37 25 14]); 
rho_max = uicontrol(param,'style','text',...
    'String','+10',...
    'Position',[230 37 25 14]); 
rho_text = uicontrol(param,'Style','text',...  
    'String','rho - distance from torus center ',...   
    'Position',[5 36 70 25]);
rho_edit = uicontrol(param,'style','edit',...
    'String',rho_slider.Value,...
    'Position',[265 34 30 25]); 
setappdata(0,'rho_slider',[rho_slider])
%% Curve Panel
curve_2D = uipanel(fig_1, 'units','normalized', 'Position',[0.75 0.20 0.25 0.6],...
    'Visible', 'On', 'Title','Intersection curve','FontSize',11);

setappdata(0,'curve_2D',[curve_2D])
setappdata(0,'fig_1',[fig_1])


ha2d = axes('Parent', curve_2D, 'Units', 'normalized', 'Position', [0.1 0.2 .8 .7]);
set(fig_1, 'CurrentAxes', ha2d)

setappdata(0,'ha2d',[ha2d])

syms x y z
f = (sqrt(x.^2 + (rho .* cosd(phi) - y .* sind(phi)).^2) - R).^2 +(z_q + y .* cosd(phi)).^2 - r.^2 

fimplicit( f == 0, 'Linewidth', 2)
daspect([1 1 1]) 

set(fig_1, 'CurrentAxes', ha3d)

%% Functions
% R slider
function R_slider_button_press(h,dummy)

    R_old = round(get(h,'Value'));
    setappdata(0,'old_R',[R_old]);
    x_old = getappdata(0,'old_x');
    y_old = getappdata(0,'old_y');
    z_old = getappdata(0,'old_z');
    design_plane(x_old, y_old, z_old)  
    
    r_old = getappdata(0,'old_r');    
    tor(R_old, r_old) 
    
    param = getappdata(0,'param');    
    R_slider = getappdata(0,'R_slider');
    R_edit = uicontrol(param,'style','edit',...
    'String',R_slider.Value,...
    'Position',[265 154 30 25]);     

    rho = getappdata(0,'old_rho');   
    phi = getappdata(0,'old_phi');
    
    fig_1 = getappdata(0,'fig_1');
    ha2d = getappdata(0,'ha2d');
    set(fig_1, 'CurrentAxes', ha2d)
    syms x y z
    f = (sqrt(x.^2 + (rho .* cosd(phi) - y .* sind(phi)).^2) - R_old).^2 +(z_old + y .* cosd(phi)).^2 - r_old.^2 
    k = fimplicit( f == 0, 'Linewidth', 2)
    XX = k.XData
    YY = k.YData
    XX = XX(~isnan(XX));
    YY = YY(~isnan(YY));
    trapz_value = trapz(XX,YY);
    
    curve_2D = getappdata(0,'curve_2D');
    r_text = uicontrol(curve_2D,'Style','text',...  % 
        'String',trapz_value ,...                % 
        'Position',[5 35 70 25]); 
    
    daspect([1 1 1]) 
    ha3d = getappdata(0,'ha3d');
    set(fig_1, 'CurrentAxes', ha3d)
end

function r_slider_button_press(h,dummy)

    r_old = round(get(h,'Value'));
    setappdata(0,'old_r',[r_old]);
    x_old = getappdata(0,'old_x');
    y_old = getappdata(0,'old_y');
    z_old = getappdata(0,'old_z');
    design_plane(x_old, y_old, z_old)  
    
    R_old = getappdata(0,'old_R');    
    tor( R_old, r_old) 

    param = getappdata(0,'param');    
    r_slider = getappdata(0,'r_slider');
    r_edit = uicontrol(param,'style','edit',...
    'String',r_slider.Value,...
    'Position',[265 124 30 25]); % L, B, W, H

    rho = getappdata(0,'old_rho');   
    phi = getappdata(0,'old_phi');
    
    fig_1 = getappdata(0,'fig_1');
    ha2d = getappdata(0,'ha2d');
    set(fig_1, 'CurrentAxes', ha2d)
    syms x y z
    f = (sqrt(x.^2 + (rho .* cosd(phi) - y .* sind(phi)).^2) - R_old).^2 +(z_old + y .* cosd(phi)).^2 - r_old.^2 
    k = fimplicit( f == 0, 'Linewidth', 2)
    XX = k.XData
    YY = k.YData
    XX = XX(~isnan(XX))
    YY = YY(~isnan(YY))
    trapz_value = trapz(XX,YY)
    
    curve_2D = getappdata(0,'curve_2D');
    r_text = uicontrol(curve_2D,'Style','text',...  % 
        'String',trapz_value ,...                % 
        'Position',[5 35 70 25]); 
    
    daspect([1 1 1]) 
    ha3d = getappdata(0,'ha3d');
    set(fig_1, 'CurrentAxes', ha3d)

end
function alpha_slider_button_press(h,dummy)
    alpha = round(get(h,'Value'));
    setappdata(0,'old_alpha',[alpha])
    
    rho = getappdata(0,'old_rho');   
    phi = getappdata(0,'old_phi');
    
    x_q = rho * cosd(alpha) * cosd(phi);
    y_q = rho * sind(alpha) * cosd(phi);
    z_q = rho * sind(phi);

    t_axis = getappdata(0,'old_t_axis')
    w_axis = getappdata(0,'old_w_axis')    
    x = x_q + t_axis * sind(alpha) - w_axis * cosd(alpha) * sind(phi);
    y = y_q - t_axis * cosd(alpha) - w_axis * sind(alpha) * sind(phi);
    z = z_q + w_axis * cosd(phi);
    
    setappdata(0,'old_x',[x]);
    setappdata(0,'old_y',[y]);
    setappdata(0,'old_z',[z]);

    design_plane(x, y, z)  
    
    R_old = getappdata(0,'old_R');   
    r_old = getappdata(0,'old_r');
    tor( R_old, r_old)
    
    fig_1 = getappdata(0,'fig_1');    
    ha2d = getappdata(0,'ha2d');   
    set(fig_1, 'CurrentAxes', ha2d)
    syms x y z
    f = (sqrt(x.^2 + (rho .* cosd(phi) - y .* sind(phi)).^2) - R_old).^2 +(z_q + y .* cosd(phi)).^2 - r_old.^2; 
    k = fimplicit( f == 0, 'Linewidth', 2)
    XX = k.XData
    YY = k.YData
    XX = XX(~isnan(XX))
    YY = YY(~isnan(YY))
    trapz_value = trapz(XX,YY)
    
    curve_2D = getappdata(0,'curve_2D');
    r_text = uicontrol(curve_2D,'Style','text',...  
        'String',trapz_value ,...                 
        'Position',[5 35 70 25]); 
    
    daspect([1 1 1]) 
    ha3d = getappdata(0,'ha3d');
    set(fig_1, 'CurrentAxes', ha3d)
    
    param = getappdata(0,'param');    
    alpha_slider = getappdata(0,'alpha_slider');
    alpha_edit = uicontrol(param,'style','edit',...
    'String', alpha_slider.Value,...
    'Position',[265 94 30 25]);    
end

function phi_slider_button_press(h,dummy)
    phi = round(get(h,'Value'));
    setappdata(0,'old_phi',[phi])
    
    rho = getappdata(0,'old_rho');   
    alpha = getappdata(0,'old_alpha');
    
    x_q = rho * cosd(alpha) * cosd(phi);
    y_q = rho * sind(alpha) * cosd(phi);
    z_q = rho * sind(phi);

    t_axis = getappdata(0,'old_t_axis')
    w_axis = getappdata(0,'old_w_axis')    
    x = x_q + t_axis * sind(alpha) - w_axis * cosd(alpha) * sind(phi);
    y = y_q - t_axis * cosd(alpha) - w_axis * sind(alpha) * sind(phi);
    z = z_q + w_axis * cosd(phi);
    
    setappdata(0,'old_x',[x]);
    setappdata(0,'old_y',[y]);
    setappdata(0,'old_z',[z]);
    
    design_plane(x, y, z)  
    
    R_old = getappdata(0,'old_R');   
    r_old = getappdata(0,'old_r');
    tor( R_old, r_old) 
    
    fig_1 = getappdata(0,'fig_1');    
    ha2d = getappdata(0,'ha2d');   
    set(fig_1, 'CurrentAxes', ha2d)
    
    syms x y z
    f = (sqrt(x.^2 + (rho .* cosd(phi) - y .* sind(phi)).^2) - R_old).^2 +(z_q + y .* cosd(phi)).^2 - r_old.^2 
    k = fimplicit( f == 0, 'Linewidth', 2)
    XX = k.XData
    YY = k.YData
    XX = XX(~isnan(XX))
    YY = YY(~isnan(YY))
    trapz_value = trapz(XX,YY)
    
    curve_2D = getappdata(0,'curve_2D');
    r_text = uicontrol(curve_2D,'Style','text',...  
        'String',trapz_value ,...                
        'Position',[5 35 70 25]); 
    
    daspect([1 1 1]) 
    ha3d = getappdata(0,'ha3d');
    set(fig_1, 'CurrentAxes', ha3d)
        
    param = getappdata(0,'param');    
    phi_slider = getappdata(0,'phi_slider');
    phi_edit = uicontrol(param,'style','edit',...
    'String',phi_slider.Value,...
    'Position',[265 64 30 25]);   

end

function rho_slider_button_press(h,dummy)    
    rho = get(h,'Value')
    setappdata(0,'old_rho',[rho]);
    
    phi = getappdata(0,'old_phi');   
    alpha = getappdata(0,'old_alpha');
    
    x_q = rho * cosd(alpha) * cosd(phi);
    y_q = rho * sind(alpha) * cosd(phi);
    z_q = rho * sind(phi);

    t_axis = getappdata(0,'old_t_axis');
    w_axis = getappdata(0,'old_w_axis');   
    x = x_q + t_axis * sind(alpha) - w_axis * cosd(alpha) * sind(phi);
    y = y_q - t_axis * cosd(alpha) - w_axis * sind(alpha) * sind(phi);
    z = z_q + w_axis * cosd(phi);
    
    setappdata(0,'old_x',[x]);
    setappdata(0,'old_y',[y]);
    setappdata(0,'old_z',[z]);
    
    design_plane(x, y, z)  
    
    R_old = getappdata(0,'old_R');   
    r_old = getappdata(0,'old_r');
    tor( R_old, r_old); 

    curve_2D = getappdata(0,'curve_2D');
    fig_1 = getappdata(0,'fig_1');
    ha2d = getappdata(0,'ha2d');
    
    set(fig_1, 'CurrentAxes', ha2d);
    syms x y z
    f = (sqrt(x.^2 + (rho .* cosd(phi) - y .* sind(phi)).^2) - R_old).^2 +(z_q + y .* cosd(phi)).^2 - r_old.^2; 
    k = fimplicit( f == 0, 'Linewidth', 2);
    XX = k.XData;
    YY = k.YData;
    XX = XX(~isnan(XX));
    YY = YY(~isnan(YY));
    trapz_value = trapz(XX,YY);

    r_text = uicontrol(curve_2D,'Style','text',...  % 
        'String',trapz_value ,...                % 
        'Position',[5 35 70 25]); 
    
    daspect([1 1 1]); 
    ha3d = getappdata(0,'ha3d');
    set(fig_1, 'CurrentAxes', ha3d);
    
    param = getappdata(0,'param');    
    rho_slider = getappdata(0,'rho_slider');
    rho_edit = uicontrol(param,'style','edit',...
    'String',rho_slider.Value,...
    'Position',[265 34 30 25]); 

end
%% call plane
function [x, y, z] = design_plane (x, y, z)

    pointA = [x(1); y(1); z(1)];
    pointB = [x(2); y(2); z(2)];
    pointC = [x(3); y(3); z(3)];

    normal = cross(pointA-pointB, pointA-pointC)
    dist = dot(normal, pointA)

    [x, y, z] = plane_surf(normal, dist, 5);
    view(50,25)
    p = surf(x, y, z);
    alpha 0.5
    hold on
end

%% create plane function helper
function [x, y, z] = plane_surf(normal, dist, size)

    normal = normal / norm(normal);
    center = normal * dist;
    tangents = null(normal') * size;

    res(1,1,:) = center + tangents * [-1;-1]; 
    res(1,2,:) = center + tangents * [-1;1]; 
    res(2,2,:) = center + tangents * [1;1]; 
    res(2,1,:) = center + tangents * [1;-1];

    x = squeeze(res(:,:,1));
    y = squeeze(res(:,:,2));
    z = squeeze(res(:,:,3));
    
end

%% create torus function

function [R, r] = tor(R,r)
    hold on
    theta_1=linspace(0,2*pi,36); % e.g. 36 partitions along perimeter of the tube 
    phi_1=linspace(0,2*pi,50); % e.g. 50 partitions along azimuth of torus

    % we convert our vectors phi and th to [n x n] matrices with meshgrid command:
    [phi_1,theta_1]=meshgrid(phi_1,theta_1); 

    x = sin(theta_1) .* (R + r*cos(phi_1));
    y = cos(theta_1) .* (R + r*cos(phi_1));
    z = r.*sin(phi_1);
    
    view(50,25)
    surface= [x; y; z];
    surf(x,y,z);
    daspect([1 1 1]);
    title('Torus')
    xlabel('X');ylabel('Y');zlabel('Z');
    alpha 0.5
    hold off
end