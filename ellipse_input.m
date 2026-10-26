%% Input function
% u = f(x_desired)
% given the semi-major and minor lengths of the ellipse. An ellipse is
% created with center as the center of the print area (dx/2, dy/2)

function [u, ref] = ellipse_input(t, a, b, T, p) % 0<=t<=T

%initializing the desired trajectory values
x = 0; y = 0;
xd = 0; yd = 0;
xdd = 0; ydd = 0;
xddd = 0; yddd = 0;

if t < T
    [q, qd, qdd, qddd] = smoothstep7(t, T);

    % ellipse parameter value
    phi = 2*pi*q;
    phid = 2*pi*qd;
    phidd = 2*pi*qdd;
    phiddd = 2*pi*qddd;

    %ellipse trajectory
    x = a*cos(phi);
    y = b*sin(phi);
    %trajectory velocity
    xd = -a*sin(phi)*phid;
    yd = b*cos(phi)*phid;
    %trajectory acceleration
    xdd = -a*(cos(phi)*phid^2 + sin(phi)*phidd);
    ydd = b*(-sin(phi)*phid^2 + cos(phi)*phidd);
    %trajectory jerk 
    xddd = a*(sin(phi)*phid^3 - 3*cos(phi)*phid*phidd - sin(phi)*phiddd);
    yddd = b*(-cos(phi)*phid^3 - 3*sin(phi)*phid*phidd + cos(phi)*phiddd);

end

% Based on the dynamics equation, currents in the left and right motors
% are calculated for the desired trajectory
i_l = (0.5*p.r/p.Km)*(p.meff*(xdd + ydd) + p.bx*xd + p.by*yd);
i_r = (0.5*p.r/p.Km)*(p.meff*(xdd - ydd) + p.bx*xd - p.by*yd);

i_ld = (0.5*p.r/p.Km)*(p.meff*(xddd + yddd) + p.bx*xdd + p.by*ydd);
i_rd = (0.5*p.r/p.Km)*(p.meff*(xddd - yddd) + p.bx*xdd - p.by*ydd);

u_left = p.Lm*i_ld + i_l*p.Rm + (p.Km/p.r)*(xd + yd);
u_right = p.Lm*i_rd + i_r*p.Rm + (p.Km/p.r)*(xd - yd);

u = [u_left, u_right]; %instantaneous input to the system

ref = struct('x', x, 'y', y ...
    , 'xdot', xd, 'ydot', yd ...
    , 'xddot', xdd, 'yddot', ydd ...
    , 'xdddot', xddd, 'ydddot', yddd);
    
end





%% smooth step function, so that the desired trajectory has no velocity, 
%acceleration and jerk, at the endpoints.
function [q, qd, qdd, qddd] = smoothstep7(t, T)
    z = min(max(t/T, 0), 1);
    q = 35*z^4 - 84*z^5 + 70*z^6 - 20*z^7;
    qd = (140*z^3 - 420*z^4 + 420*z^5 - 140*z^6)/T;
    qdd = (420*z^2 - 1680*z^3+ 2100*z^4 - 840*z^5)/T^2;
    qddd = (840*z - 5040*z^2 + 8400*z^3 - 4200*z^4)/T^3;

end