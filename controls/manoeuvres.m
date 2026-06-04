function control = manoeuvres(trajectory, sim)
%MANOEUVRES Define waypoint trajectories for the XDrone simulation.
%
% Options:
%   'figure8' | 'hstep' | 'vstep' | 'straight' | 'down'

    R  = 30;                  % circle radius [m]
    z0 = sim.initial.pos0(3); % altitude [m, NED]

    loop = true;

    switch trajectory

        case 'figure8'
            wp = [
                R*cosd(315),  R + R*sind(315),  z0;
                R*cosd(  0),  R + R*sind(  0),  z0;
                R*cosd( 45),  R + R*sind( 45),  z0;
                R*cosd( 90),  R + R*sind( 90),  z0;
                R*cosd(135),  R + R*sind(135),  z0;
                R*cosd(180),  R + R*sind(180),  z0;

                R*cosd(  0), -R + R*sind(  0),  z0;
                R*cosd(315), -R + R*sind(315),  z0;
                R*cosd(270), -R + R*sind(270),  z0;
                R*cosd(225), -R + R*sind(225),  z0;
                R*cosd(180), -R + R*sind(180),  z0;
                R*cosd(135), -R + R*sind(135),  z0;
            ];

        case 'hstep'
            wp = [
                 50,    0,  z0;
                 50,  150,  z0;
                200,  150,  z0;
                260,  150,  z0;
            ];
            loop = false;

        case 'up'
            wp = [
                 50,   0,  z0;
                 50,   0,  z0-150;
                200,   0,  z0-150;
                
            ];
            loop = false;

        case 'down'
            wp = [
                50,   0,  z0;
                300,  0,  z0+50;
                400,  0,  z0+50;
                
            ];
            loop = false;

        case 'straight'
            wp = [
                100,  0,  z0;
                150,  0,  z0;
            ];
            loop = false;

        otherwise
            error('Unknown trajectory ''%s''. Choose: figure8 | hstep | vstep | straight | down', trajectory);
    end

    control.waypoints = wp';
    control.waypoint_radius = 30.0;
    control.loop = loop;
end