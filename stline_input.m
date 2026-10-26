function [u, ref] = stline_input(x1, y1, x2, y2, T, t, p)
% Feedforward input for a straight-line move from (x1,y1) to (x2,y2)
x = x2; y = y2;             % hold at the endpoint once finished
xd = 0; yd = 0;
xdd = 0; ydd = 0;
xddd = 0; yddd = 0;

if t < T
    [q, qd, qdd, qddd] = smoothstep7(t, T);
    x    = x1*(1-q) + x2*q;
    y    = y1*(1-q) + y2*q;
    xd   = (x2 - x1)*qd;    yd   = (y2 - y1)*qd;
    xdd  = (x2 - x1)*qdd;   ydd  = (y2 - y1)*qdd;
    xddd = (x2 - x1)*qddd;  yddd = (y2 - y1)*qddd;
end

% Motor currents and their derivatives from the mechanical dynamics
i_l  = (0.5*p.r/p.Km)*(p.meff*(xdd  + ydd)  + p.bx*xd  + p.by*yd);
i_r  = (0.5*p.r/p.Km)*(p.meff*(xdd  - ydd)  + p.bx*xd  - p.by*yd);
i_ld = (0.5*p.r/p.Km)*(p.meff*(xddd + yddd) + p.bx*xdd + p.by*ydd);
i_rd = (0.5*p.r/p.Km)*(p.meff*(xddd - yddd) + p.bx*xdd - p.by*ydd);

% Motor voltages
u_left  = p.Lm*i_ld + p.Rm*i_l + (p.Km/p.r)*(xd + yd);
u_right = p.Lm*i_rd + p.Rm*i_r + (p.Km/p.r)*(xd - yd);
u = [u_left, u_right];

ref = struct('x', x, 'y', y, ...
             'xdot', xd, 'ydot', yd, ...
             'xddot', xdd, 'yddot', ydd, ...
             'xdddot', xddd, 'ydddot', yddd);
end

function [q, qd, qdd, qddd] = smoothstep7(t, T)
% 7th-order smoothstep: zero vel/acc/jerk at both endpoints
    z = min(max(t/T, 0), 1);
    q    = 35*z^4 - 84*z^5 + 70*z^6 - 20*z^7;
    qd   = (140*z^3 - 420*z^4 + 420*z^5 - 140*z^6)/T;
    qdd  = (420*z^2 - 1680*z^3 + 2100*z^4 - 840*z^5)/T^2;
    qddd = (840*z - 5040*z^2 + 8400*z^3 - 4200*z^4)/T^3;
end