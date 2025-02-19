function [t,y] = solve_slow_ODE(G,B,V_0,V_1,delta,J0_1A,J1_1A,J0_2A,t_span,Init,options)
%SOLVE_SLOW_ODE Solve the emergent system of 6 ODEs governing the rotation and
%translation of a rapidly yawing spheroid.


[t,y] = ode15s(@(t,y) ode_num_slow(t,y,G,B,V_0,V_1,delta,J0_1A,J1_1A,J0_2A),t_span,Init,options);

    function dy = ode_num_slow(t,y,G,B,V_0,V_1,delta,J0_1A,J1_1A,J0_2A)

        dy = zeros(size(y));
        
        alp = y(1);
        bet = y(2);
        gam = y(3);
        
        f1_B = -B*sin(alp)*((sin(bet)^2 + J0_2A*(1 + cos(bet)^2))*cos(alp)*sin(2*gam) ...
            - (1 - J0_2A)*sin(bet)*cos(bet)*cos(2*gam));
        f2_B = B*((1 - J0_2A)*sin(bet)*cos(bet)*sin(2*gam) ...
            + (sin(bet)^2 + J0_2A*(1 + cos(bet)^2))*cos(alp)*cos(2*gam));
        f3_B = 1 - (B/2)*((cos(bet)^2 + J0_2A*(1 + sin(bet)^2))*cos(2*gam) ...
            - (1 - J0_2A)*sin(bet)*cos(alp)*cos(bet)*sin(2*gam));

        % Rotational dynamics
        dy(1) = G*(f1_B/4);
        dy(2) = G*(f2_B/4);
        dy(3) = G*(f3_B/2);

        e1hat_vector = [cos(alp), sin(alp)*sin(gam), -sin(alp)*cos(gam)];
        e2hat_vector = [sin(alp)*sin(bet), cos(gam)*cos(bet) - sin(gam)*cos(alp)*sin(bet), ...
            sin(gam)*cos(bet) + cos(gam)*cos(alp)*sin(bet)];
        e3hat_vector = [sin(alp)*cos(bet), -cos(gam)*sin(bet) - sin(gam)*cos(alp)*cos(bet), ...
            -sin(gam)*sin(bet) + cos(gam)*cos(alp)*cos(bet)];

        % Translational dynamics
        dy(4) = (V_0(1)*J0_1A + V_1(3)*J1_1A*sin(delta(3)))*e1hat_vector(1) ...
            + (V_0(3)*J0_1A - V_1(1)*J1_1A*sin(delta(1)))*e3hat_vector(1) ...
            + V_0(2)*e2hat_vector(1);
        dy(5) = (V_0(1)*J0_1A + V_1(3)*J1_1A*sin(delta(3)))*e1hat_vector(2) ...
            + (V_0(3)*J0_1A - V_1(1)*J1_1A*sin(delta(1)))*e3hat_vector(2) ...
            + V_0(2)*e2hat_vector(2);
        dy(6) = (V_0(1)*J0_1A + V_1(3)*J1_1A*sin(delta(3)))*e1hat_vector(3) ...
            + (V_0(3)*J0_1A - V_1(1)*J1_1A*sin(delta(1)))*e3hat_vector(3) ...
            + V_0(2)*e2hat_vector(3) + G*y(5);

    end

end

