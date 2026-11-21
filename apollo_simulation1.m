clear; clc; close all;

% constants
G  = 6.67430e-11;  % gravitational constant (SI units)

% masses of the earth and moon in kg
m_e = 5.9722e24;  % mass of the earth
m_m = 7.34203e22;  % mass of the moon

% spacecraft mass before and after the S-IVB separation (kg)
m_stack0 = 140121.9;  % mass with S-IVB attached
m_sc0    = 45182;  % mass after dropping S-IVB

% constants for calculating mass lost during the burns
g0 = 9.80665;  % standard gravity (m/s^2)
IspTLI = 421;  % specific impulse for the TLI burn (s)
IspSPS = 314;  % specific impulse for the later burns (s)
veTLI = IspTLI*g0;  % exhaust velocity (m/s)
veSPS = IspSPS*g0;  % exhaust velocity (m/s)

% set the timestep and simulation length
dt = 5;  % timestep in seconds
tmax = 3.4*24*3600+6*3600;  % maximum simulation time in seconds
clockmax = floor(tmax/dt);

% initial position and velocity of the moon from JPL Horizons
% the original coordinates are in km, so convert to meters
x_m0  = -2.892936360244735E+05 * 1e3;
y_m0  =  2.810455841981057E+05 * 1e3;
vx_m0 = -6.654083626094236E-01 * 1e3;
vy_m0 = -7.110442245926685E-01 * 1e3;

% set the spacecraft starting position and velocity in parking orbit
rpark = 6556e3;  % distance from the center of earth (m)
vpark = -sqrt((G*m_e)/rpark);  % parking orbit velocity (m/s)

theta0 = 0;  % starting angle (radians)
x_sc0  = rpark*cos(theta0);
y_sc0  = rpark*sin(theta0);
vx_sc0 = -vpark*sin(theta0);
vy_sc0 =  vpark*cos(theta0);

% times, magnitudes, and angles of the scheduled burns
t_TLI = (2*3600 + 44*60 + 16.2 +650);  % TLI time in seconds
t_SEP = t_TLI + 3600;  % separate the S-IVB one hour after TLI
t_MCC = (26*3600 + 44*60 + 58.64);  % midcourse correction time (s)

dV_TLI = 3264;  % TLI change in speed (m/s)
dV_MCC = 50;  % MCC change in speed (m/s)

TLIangle = -18.114;  % TLI angle in degrees
MCCangle = 10;  % MCC angle in degrees

% conditions for entering and circularizing lunar orbit
Rmoon = 1.7374e6;  % radius of the moon (m)
rp_min = Rmoon + 110e3;  % minimum radius for the lunar orbit
r_capture_trigger = 9.0e7;  % distance where the orbit insertion check begins (m)
t_wait_circ = 4*3600;  % time between LOI and circularization
extra_time_after_circ = 8*3600;  % extra time after circularization before stopping

% store the masses, positions, and velocities of the three bodies
% body 1 is Earth, 2 is the moon, and 3 is the spacecraft
n = 3;

m = [m_e; m_m; m_stack0];  % masses of the three bodies (kg)
x = [0; x_m0; x_sc0];  % x positions (m)
y = [0; y_m0; y_sc0];  % y positions (m)
u = [0; vx_m0; vx_sc0];  % x velocities (m/s)
v = [0; vy_m0; vy_sc0];  % y velocities (m/s)

% arrays to store values during the simulation
tsave  = zeros(1,clockmax);
xsave  = zeros(n,clockmax);
ysave  = zeros(n,clockmax);
msave  = zeros(1,clockmax);
Esave  = zeros(1,clockmax);  % store total energy
dEsave = zeros(1,clockmax);  % store change in energy
rMsave = zeros(1,clockmax);

% set up the trajectory plot
cEarth=[0.45 0.75 1.0];  % earth is blue
cMoon =[1 1 1];  % moon is white
cSC   =[1.0 0.35 0.75];  % spacecraft is pink
cBurn =[1.0 0.55 0.0];  % burn markers are orange

figure(1); clf; hold on; axis equal; set(gca,'Color','k');

a = 1.5*hypot(x_m0,y_m0);  % plotting range
axis([-a a -a a]);

hp_e = plot(x(1),y(1),'o','Color',cEarth,'MarkerFaceColor',cEarth,'MarkerSize',10,'LineWidth',2);
hp_m = plot(x(2),y(2),'o','Color',cMoon ,'MarkerFaceColor',cMoon ,'MarkerSize',10,'LineWidth',2);
hp_s = plot(x(3),y(3),'o','Color',cSC   ,'MarkerFaceColor',cSC   ,'MarkerSize',8 ,'LineWidth',1.5);

ht_m = plot(x(2),y(2),'-','Color',cMoon,'LineWidth',1.5);
ht_s = plot(x(3),y(3),'-','Color',cSC  ,'LineWidth',1.8);

title('Trajectory (Earth–Moon–Spacecraft)','Color','w');
xlabel('x (m)','Color','w'); ylabel('y (m)','Color','w');
set(gca,'XColor','w','YColor','w');

% set up the animation recording
makeVideo = true;  % set to false if no video is needed
vidName   = 'apollo_trajectory_project.mp4';  % name of the animation file
vidFPS    = 30;  % frames per second in the animation

if makeVideo
    set(gcf,'Color','k');  % black background for the video
    vobj = VideoWriter(vidName,'MPEG-4');
    vobj.FrameRate = vidFPS;
    vobj.Quality   = 100;  % highest video quality
    open(vobj);
end

% boolean variables to make sure each burn only happens once
didTLI=false; didSEP=false; didMCC=false;
entered_capture=false; didLOI=false; didCIRC=false;
t_CIRC_target = NaN;

% previous moon distances to check for the closest approach
rM_prev2 = NaN;
rM_prev1 = NaN;

% calculate the initial total energy of the system
E0 = 0;

% kinetic energy
for i=1:n
    E0 = E0 + 0.5*m(i)*(u(i)^2 + v(i)^2);
end

% potential energy
for i=1:n
    for j=i+1:n
        rij = hypot(x(j)-x(i), y(j)-y(i));
        E0 = E0 - G*m(i)*m(j)/rij;
    end
end

% main loop using the Symplectic Euler method
for clock = 1:clockmax

    t = (clock-1)*dt;

    % store values
    tsave(clock) = t;
    xsave(:,clock) = x;
    ysave(:,clock) = y;
    msave(clock) = m(3);

    % calculating total energy
    E = 0;

    % kinetic energy
    for i=1:n
        E = E + 0.5*m(i)*(u(i)^2 + v(i)^2);
    end

    % potential energy
    for i=1:n
        for j=i+1:n
            rij = hypot(x(j)-x(i), y(j)-y(i));
            E = E - G*m(i)*m(j)/rij;
        end
    end

    Esave(clock)  = E;
    dEsave(clock) = E - E0;

    % spacecraft distance to Moon
    rM = hypot(x(3)-x(2), y(3)-y(2));
    rMsave(clock) = rM;

    % for loop to update velocities via Symplectic Euler
    for i = 1:n
        for j = 1:n
            if i ~= j
                dx = x(j) - x(i);
                dy = y(j) - y(i);
                r  = sqrt(dx*dx + dy*dy);

                u(i) = u(i) + dt * G*m(j)*dx/(r^3);
                v(i) = v(i) + dt * G*m(j)*dy/(r^3);
            end
        end
    end

    %Burns section
    %store the current spacecraft and moon velocity for easier delta v computations
    v_sc = [u(3); v(3)];
    v_mo = [u(2); v(2)];
    vrelM = v_sc - v_mo;

    % TLI Burn
    if ~didTLI && t >= t_TLI

        %direction is unit vector in direction of velocity
        vdir = v_sc / norm(v_sc);
        %convert the burn angle to radians for the rotation matrix
        burnangle = deg2rad(TLIangle);
        %2D rotation matrix to aim the burn
        R = [cos(burnangle) -sin(burnangle); sin(burnangle) cos(burnangle)];
        %compute direction by rotating vdir
        dir = R * vdir;

        %burn magnitude set
        dv = dV_TLI;

        %calculate mass ratio from rocket equation
        mratio = exp(-dv/veTLI);
        m(3) = m(3) * mratio;

        %change velocity components by burn delta v
        u(3) = u(3) + dv*dir(1);
        v(3) = v(3) + dv*dir(2);

        didTLI = true;

        % plot the burn location
        plot(x(3),y(3),'o','Color',cBurn,'MarkerSize',14,'LineWidth',2.8);
    end

    % Separation of S-IVB
    if ~didSEP && t >= t_SEP
        %drop mass by the mass of the S-IVB
        m(3) = m_sc0;
        didSEP = true;

        %no burn marker for the S-IVB separation
        %plot(x(3),y(3),'o','Color',cBurn,'MarkerSize',14,'LineWidth',2.8);
    end

    % Midcourse correction
    if ~didMCC && t >= t_MCC

        %direction is unit vector in direction of velocity
        vdir = v_sc / norm(v_sc);
        %convert the burn angle to radians for the rotation matrix
        burnangle = deg2rad(MCCangle);
        %2D rotation matrix to aim the burn
        R = [cos(burnangle) -sin(burnangle); sin(burnangle) cos(burnangle)];
        %compute direction by rotating vdir
        dir = R * vdir;

        %burn magnitude set
        dv = dV_MCC;

        %calculate mass ratio from rocket equation
        mratio = exp(-dv/veSPS);
        m(3) = m(3) * mratio;

        %change velocity components by burn delta v
        u(3) = u(3) + dv*dir(1);
        v(3) = v(3) + dv*dir(2);

        didMCC = true;

        % plot the burn location
        plot(x(3),y(3),'o','Color',cBurn,'MarkerSize',14,'LineWidth',2.8);
    end

    % detects when the spacecraft is close enough to the moon for orbit insertion
    if ~entered_capture && rM < r_capture_trigger
        entered_capture = true;
    end

    % Lunar orbit insertion, once the spacecraft reaches its closest approach
    if entered_capture && ~didLOI

        %make sure there is enough history to check for a minimum distance
        if ~isnan(rM_prev2) && ~isnan(rM_prev1)

            %the LOI burn takes place at a local min in distance to the moon
            if (rM_prev1 < rM_prev2) && (rM_prev1 < rM)

                %set burn radius, but make sure it is not less than the minimum
                rloi = max(rM_prev1, rp_min);

                vrelM = [u(3)-u(2); v(3)-v(2)];
                normvrelM = norm(vrelM);

                ra = rloi + 203720;  % difference between desired apolune and perilune radii
                %construct the desired elliptical orbit around the moon
                ades = 0.5*(rloi + ra);

                vdes = sqrt((G*m_m)*(2/rloi - 1/ades));

                %burn magnitude, never allow this burn to speed the spacecraft up
                dv = max(0, normvrelM - vdes);

                if dv > 0
                    %burn opposite the velocity relative to the moon to slow down
                    dir = -vrelM / norm(vrelM);

                    %calculate mass ratio from rocket equation
                    mratio = exp(-dv/veSPS);
                    m(3) = m(3) * mratio;

                    %change velocity components by burn delta v
                    u(3) = u(3) + dv*dir(1);
                    v(3) = v(3) + dv*dir(2);
                end

                didLOI = true;
                t_CIRC_target = t + t_wait_circ;

                % plot the burn location
                plot(x(3),y(3),'o','Color',cBurn,'MarkerSize',14,'LineWidth',2.8);
            end
        end

        rM_prev2 = rM_prev1;
        rM_prev1 = rM;
    end

    % circularization burn
    if didLOI && ~didCIRC && t >= t_CIRC_target

        r = hypot(x(3)-x(2), y(3)-y(2));
        vrelM = [u(3)-u(2); v(3)-v(2)];
        normvrelM = norm(vrelM);

        %needed circular orbit speed
        vcirc = sqrt((G*m_m)/r);
        %required burn magnitude
        dv = abs(vcirc - normvrelM);

        if dv > 0
            %check if the spacecraft is too slow or too fast and burn accordingly
            if vcirc > normvrelM
                dir = vrelM / norm(vrelM);
            else
                %burn opposite the velocity relative to the moon to slow down
                dir = -vrelM / norm(vrelM);
            end

            %calculate mass ratio from rocket equation
            mratio = exp(-dv/veSPS);
            m(3) = m(3) * mratio;

            %change velocity components by burn delta v
            u(3) = u(3) + dv*dir(1);
            v(3) = v(3) + dv*dir(2);
        end

        didCIRC = true;

        % plot the burn location
        plot(x(3),y(3),'o','Color',cBurn,'MarkerSize',14,'LineWidth',2.8);
    end

    % Updating position values, complete
    for i = 1:n
        x(i) = x(i) + dt*u(i);
        y(i) = y(i) + dt*v(i);
    end

    % update plots and write video frames at the selected timesteps
    if mod(clock,120)==0

        hp_e.XData = x(1); hp_e.YData = y(1);
        hp_m.XData = x(2); hp_m.YData = y(2);
        hp_s.XData = x(3); hp_s.YData = y(3);

        ht_m.XData = xsave(2,1:clock); ht_m.YData = ysave(2,1:clock);
        ht_s.XData = xsave(3,1:clock); ht_s.YData = ysave(3,1:clock);

        drawnow;

        % save the current trajectory plot to the video
        if makeVideo
            frame = getframe(gcf);
            writeVideo(vobj, frame);
        end
    end

    %stop a short time after circularization
    if didCIRC && t > (t_CIRC_target + extra_time_after_circ)
        break;
    end
end

% close the video after the simulation is finished
if makeVideo
    close(vobj);
    fprintf("Saved MP4 animation: %s\n", vidName);
end

% remove the unused end of the storage arrays
kmax = clock;
tsave  = tsave(1:kmax);
msave  = msave(1:kmax);
Esave  = Esave(1:kmax);
dEsave = dEsave(1:kmax);

fprintf("Flags: TLI=%d SEP=%d MCC=%d LOI=%d CIRC=%d\n",didTLI,didSEP,didMCC,didLOI,didCIRC);
fprintf("Closest approach achieved: %.0f km\n", min(rMsave(1:kmax))/1000);

% plot total energy over time
figure(2); clf;
plot(tsave, Esave, 'LineWidth', 1.6); grid on;
xlabel('Time (s)');
ylabel('Total Energy E(t) (J)');
title('Total Energy vs Time');

Emean = mean(Esave);
window = abs(Emean) * 1e-3;  % set the plotted energy range
ylim([Emean - window, Emean + window]);

% plot relative change in total energy
figure(4); clf;

% relative energy change compared to the initial total energy
dErel = dEsave ./ abs(E0);

plot(tsave, dErel, 'LineWidth', 1.6);
grid on;

xlabel('Time (s)');
ylabel('\DeltaE / |E_0|');
title('Relative Energy Drift');

% keep the y axis centered around zero
dr = max(abs(dErel));

if dr == 0
    dr = 1e-16;
end

ylim([-1.1*dr, 1.1*dr]);

% show the y axis in scientific notation
ax = gca;
ax.YAxis.Exponent = floor(log10(dr));

% plot spacecraft mass over time
figure(3); clf;
plot(tsave, msave, 'LineWidth', 1.8); grid on;
xlabel('Time (s)'); ylabel('Vehicle mass (kg)');
title('Vehicle Mass vs Time');

mmin=min(msave); mmax=max(msave);
ylim([mmin-0.02*(mmax-mmin), mmax+0.02*(mmax-mmin)]);