%% Input function for a Lissajous curve centered at the origin
% x = a*sin(nx*phi + delta),  y = b*sin(ny*phi),  phi = 2*pi*q(t)
% nx, ny: integer frequency ratio, delta: phase shift (pi/2 starts at (a,0))
function [u, ref] = lissajous_input(t, a, b, nx, ny, delta, T, p) % 0<=t<=T
    [q, qd, qdd, qddd] = smoothstep7(t, T);   

    % curve parameter and its derivatives
    phi    = 2*pi*q;
    phid   = 2*pi*qd;
    phidd  = 2*pi*qdd;
    phiddd = 2*pi*qddd;

    th_x = nx*phi + delta;
    th_y = ny*phi;

    % trajectory
    x = a*sin(th_x);
    y = b*sin(th_y);
    % velocity
    xd = a*nx*cos(th_x)*phid;
    yd = b*ny*cos(th_y)*phid;
    % acceleration
    xdd = a*nx*(-nx*sin(th_x)*phid^2 + cos(th_x)*phidd);
    ydd = b*ny*(-ny*sin(th_y)*phid^2 + cos(th_y)*phidd);
    % jerk
    xddd = a*nx*(-nx^2*cos(th_x)*phid^3 - 3*nx*sin(th_x)*phid*phidd + cos(th_x)*phiddd);
    yddd = b*ny*(-ny^2*cos(th_y)*phid^3 - 3*ny*sin(th_y)*phid*phidd + cos(th_y)*phiddd);

    % current from dynamics
    i_l  = (0.5*p.r/p.Km)*(p.meff*(xdd  + ydd)  + p.bx*xd  + p.by*yd);
    i_r  = (0.5*p.r/p.Km)*(p.meff*(xdd  - ydd)  + p.bx*xd  - p.by*yd);
    i_ld = (0.5*p.r/p.Km)*(p.meff*(xddd + yddd) + p.bx*xdd + p.by*ydd);
    i_rd = (0.5*p.r/p.Km)*(p.meff*(xddd - yddd) + p.bx*xdd - p.by*ydd);

    % voltages 
    u_left  = p.Lm*i_ld + p.Rm*i_l + (p.Km/p.r)*(xd + yd);
    u_right = p.Lm*i_rd + p.Rm*i_r + (p.Km/p.r)*(xd - yd);
    u = [u_left, u_right];

    ref = struct('x', x, 'y', y, ...
                 'xdot', xd, 'ydot', yd, ...
                 'xddot', xdd, 'yddot', ydd, ...
                 'xdddot', xddd, 'ydddot', yddd);
end

% smoothstep7 fxn for smoother trajectory
function [q, qd, qdd, qddd] = smoothstep7(t, T)
    z = min(max(t/T, 0), 1);
    q    = 35*z^4 - 84*z^5 + 70*z^6 - 20*z^7;
    qd   = (140*z^3 - 420*z^4 + 420*z^5 - 140*z^6)/T;
    qdd  = (420*z^2 - 1680*z^3 + 2100*z^4 - 840*z^5)/T^2;
    qddd = (840*z - 5040*z^2 + 8400*z^3 - 4200*z^4)/T^3;
end
