%%% Austin Kilgore - Aerodynamics Project - NACA Thin Airfoil

clear
clc



%%% User Inputs

naca = input("Enter 4-Digit NACA Number: ",'s'); % the four-digit NACA designation; 
% please enter as a character array (use single quotes)

aoa = input("Enter Angle of Attack (in degrees): "); 
% angle of attack; enter value in degrees; the limit is ~10

fprintf("Note: Process Takes ~10 seconds (depending on parfor and if" + ...
    " program was used recently)\n") % info for user



%%% Part I: Determine Lift Information

tic % time the process

% process input data

max_cam = str2double(naca(1)) / 100; % max camber

peak_loc = str2double(naca(2)) / 10; % location of max camber

peak_loc_theta = acos(1 - 2*peak_loc); % location of max camber after transform

max_thick = str2double(naca(3:4)) / 100; % maximum thickness

aoa = deg2rad(aoa); % switch to radians

% computing zero-lift angle of attack

l0aoa = (-1/pi)*(((2*max_cam/(peak_loc^2))*(peak_loc_theta*(0.75 - peak_loc) +...
        sin(peak_loc_theta)*(peak_loc + 0.25*cos(peak_loc_theta) - 1)))  + ...
        (2*max_cam/((1-peak_loc)^2)*(peak_loc*(peak_loc_theta - pi) - ...
        sin(peak_loc_theta)*(peak_loc + 0.25*cos(peak_loc_theta) - 1) - ...
        0.75*peak_loc_theta + 0.75*pi))); % painful to write

% credit to Wolfram for evaluating the two piecewise integrals (and others)

if max_cam == 0
    l0aoa = 0; % symmetric failsafe, to avoid NaN
end

% sectional lift coefficient

cl = 2*pi*(aoa - l0aoa); % basic thin airfoil lift equation   



%%% Part II: Create Airfoil Visual
% note: all further parts of this project are for the visual representation 
% of the flow; the sectional lift coefficient and the zero-lift angle of  
% attack has already been determined

% thickness parameters and general equations from airfoiltools.com

a0 = 0.2969;
a1 = -.126;
a2 = -0.3516;
a3 = .2843;
a4 = -.1036;

% initializing arrays and parameters

space = 0.005; % mesh size

divpt = peak_loc/space; % the piecewise break (index)

dist = 0:space:1; % effectively the 'x' vector

thick = (max_thick/.2)*(a0*(dist.^(0.5)) + a1.*dist + a2*(dist.^2) + ...
            a3*(dist.^3) + a4*(dist.^4)); % half-thickness along the chord,
% again, credit to airfoil tools for the specific equations

camber = zeros(size(dist)); % camber line
suct = zeros(size(dist)); % suction side
pres = zeros(size(dist)); % pressure side

% calculating airfoil values (camber, suction, and pressure)
% essentially, creating a camber vector, then adding or 
% subtracting the half-thickness from above to get suction and pressure

for i = 1:numel(dist)
    if max_cam ~= 0 % symmetric failsafe
        if i <= (divpt + 1)
            camber(i) = (max_cam/(peak_loc^2))*(2*peak_loc*dist(i) - dist(i)^2);
            % leading portion of piecewise equation
        else
            camber(i) = (max_cam/((1 - peak_loc)^2))*...
                (1 - 2*peak_loc + 2*peak_loc*dist(i) - dist(i)^2);
            % trailing edge of "                  "
        end
    end
    suct(i) = camber(i) - thick(i); % just use camber and half-thickness
    pres(i) = camber(i) + thick(i); % to calculate suction and pressure
end



%%% Part III: Potential and Velocity

% x and y limits of plot (which leads to calculations, so don't go insane)

xl = [-2 3]; 
yl = [-2 2]; 

% fourier coefficients 

flim = 500; % number of desired fourier coefficients plus one 
% (for the circulation per unit length equation)

if max_cam == 0
    a_0 = aoa; % symmetric failsafe
else
a_0 = aoa - (1/pi)*(((2*max_cam/(peak_loc^2))*((peak_loc-0.5)*...
    peak_loc_theta+0.5*sin(peak_loc_theta)))+...
    ((2*max_cam/((1-peak_loc)^2))*(peak_loc*(pi-peak_loc_theta)+...
    0.5*peak_loc_theta-0.5*sin(peak_loc_theta)-(pi/2)))); % A_0
end

a_n = linspace(0,0,flim); % initialize vector for A_n's

for n = 1:flim % A_n's, n ~= 0
    if n == 1 % due to singularities in general formula, separate A_1
        a_n(n) = (2/pi)*(((2*max_cam/(peak_loc^2))*(sin(peak_loc_theta)*...
            (peak_loc+0.25*cos(peak_loc_theta)-0.5)+0.25*peak_loc_theta)) + ...
            ((2*max_cam/((1-peak_loc)^2))*(-sin(peak_loc_theta)*...
            (peak_loc+0.25*cos(peak_loc_theta)-0.5)-0.25*peak_loc_theta+0.25*pi)));
    else % other A_n's
        a_n(n) = (2/pi)*(((2*max_cam/(peak_loc^2))*...
            ((((sin(n*peak_loc_theta))*(((0.5*(n^2)*...
            cos(peak_loc_theta))/((n^2)-1))+peak_loc-0.5))/n)-...
            ((0.5*sin(peak_loc_theta)*cos(n*peak_loc_theta))/((n^2)-1))))+...
            ((((2*max_cam)/((1-peak_loc)^2))*(1/n/((n^2)-1)))*...
            ((sin(n*peak_loc_theta)*(((n^2)*(0.5-peak_loc))-...
            (0.5*(n^2)*cos(peak_loc_theta))+peak_loc-0.5))+...
            ((n^2)*(peak_loc-1)-peak_loc+0.5)+(0.5*n*sin(peak_loc_theta)*...
            cos(n*peak_loc_theta))))); % credit to wolfram for integrals
    end
    if isnan(a_n(n))
        a_n(n) = 0; % symmetric failsafe
    end
end

% vortex calculations

num_vort = 2000; % number of vortices for GRAPHICAL representation

dtheta = pi/(num_vort + 1); % spacing of vortices in theta-space

dgamma = linspace(0, 0, num_vort); % vector of gammas at each vortex 
% (discretization of circulation per unit length function)

thetapos = linspace(0.5*dtheta, pi-(0.5*dtheta), num_vort); 
% vortex position in theta-space, I purposely avoided the boundaries of the
% interval since I expect that issues would arise

for n = 1:num_vort % calculating the dgamma's
    dg = 0;
    for i = 1:flim % the inner part of the double sum (A_n's)
        dg = dg + a_n(i)*sin(i*thetapos(n));
    end
    dg = dg*sin(thetapos(n)); % finishing up the expression
    dg = dg + a_0*(1 + cos(thetapos(n))); 
    dg = dg*dtheta*(1/(2*pi)); 
    dgamma(n) = dg; 
    % strictly speaking, this isn't gamma since I included the /2/pi, etc
    % the purpose is so the stream function is this times ln(r)
    % I use the stream function instead of the potential function due to
    % the discontinuity at pi for arctan2; believe me, I tried
end

% mesh prep

xpos = 0.5*(1-cos(thetapos)); % positions of vortices in x-space

x = xl(1):space:xl(2);
y = yl(1):space:yl(2); % using previously stated x and y limits

xsteps = numel(x);
ysteps = numel(y); % useful later

[x,y] = meshgrid(x,y); % need a matrix of x and y values

% stream function calculations

stream = (cos(aoa).*y - sin(aoa).*x); % free-stream

parfor n = 1:num_vort % stream func induced from each mini-vortex
    stream = stream + dgamma(n)*log(sqrt(y.^2 + (x-xpos(n)).^2));
end
% NOTE: The above loop takes up nearly all of the runtime of this program
% So, if you need to reduce it significantly, change num_vort
% This is also why parfor was used

% velocity calculations 

dPsidx = zeros(ysteps, xsteps);
dPsidy = zeros(ysteps, xsteps);

for i = 2:(ysteps-1) % taken from my HW 5, finite difference
    for j = 2:(xsteps - 1) 
        dPsidx(i,j) = (stream(i,j+1) - stream(i,j-1)) / (2*space);
        dPsidy(i,j) = (stream(i+1,j) - stream(i-1,j)) / (2*space);
    end
end

u = dPsidy; % getting velocities from stream funcs
v = -dPsidx;

u_tot = sqrt(u.^2 + v.^2); % total velocity

% using NaN to get rid of velocity values inside of the airfoil

x0 = ((-xl(1)/space)+1); % indices of x=0 and x=1 (LE and TE)
x1 = (((1-xl(1))/space)+1);

for i = x0:x1 % find indices in u_tot that are inside airfoil, erase slice
    li = (((yl(2)-yl(1))/(2*space)) + 1 + round((pres(i+(xl(1)/space))/space)));
    ui = (((yl(2)-yl(1))/(2*space)) + 1 + round((suct(i+(xl(1)/space))/space)));
    u_tot(ui:li,i) = NaN;
    u(ui:li,i) = NaN;
    v(ui:li,i) = NaN;
end



%%% Part IV: Plotting

% speed/stream 3D plot

daspect([1 1 1]) % prevent dialation

vis00 = surf(u_tot); % gives a color/elevation visual of the speed of the flow
shading interp;

vis0 = streamslice(u,v,5,'noarrows'); % streamlines

set(vis0,'color','k');

for i=1:length(vis0)
    zi = interp2(u_tot,vis0(i).XData, vis0(i).YData); 
    vis0(i).ZData = zi;
    % this places the streamlines onto the speed plot, pretty cool
end

colorbar
caxis([0.5 2]) %#ok<CAXIS>
zlim([0.5 2]) % takes out bad values
xlabel('Axis Parallel to Free-Stream')
ylabel('Axis Normal to Free-Stream')
zlabel('Speed of Flow')
rotate(vis00, [0 0 1], -rad2deg(aoa), [0 0 0]) % rotate so free-stream
rotate(vis0, [0 0 1], -rad2deg(aoa), [0 0 0]) % is parallel to an axis

title('Local Speed of Air Around the Airfoil') 

% standard 2D stream plot with airfoil

figure

hold on

xlim(xl) % limits of graph
ylim(yl)

vis2 = streamslice(x,y,u,v,5,'noarrows'); % streamlines

vis1 = plot(dist, suct, dist, pres); % outline of airfoil

x2 = [dist, fliplr(dist)]; % credit for filling code to 'Image Analyst'
inBetween = [pres, fliplr(suct)]; % https://www.mathworks.com/matlabcentral
fill = fill(x2, inBetween, 'g'); % /answers/180829-shade-area-between-graphs

rotate(vis1, [0 0 1], -rad2deg(aoa), [0 0 0]) % rotate the figure such that
rotate(vis2, [0 0 1], -rad2deg(aoa), [0 0 0]) % it's in a 'normal' frame
rotate(fill, [0 0 1], -rad2deg(aoa), [0 0 0])

daspect([1 1 1]) % prevent dialation

title("Zero Lift Angle of Attack: " + num2str(rad2deg(l0aoa)) + char(176) + ...
    newline + "Sectional Lift Coefficient: " + num2str(cl)) % present data

hold off

toc

