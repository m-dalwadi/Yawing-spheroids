function [t,y] = solve_ODE(G,A,B,V_0,V_1,delta,w,T_init,T_max,Init,options)
%SOLVE_ODE Solve the full system of 6 ODEs governing the rotation and
%translation of a rapidly yawing spheroid.


[t,y] = ode15s(@(t,y) ode_num(t,y,G,A,B,V_0,V_1,delta,w),[T_init T_max],Init,options);


    function dy = ode_num(t,y,G,A,B,V_0,V_1,delta,w)

        dy = zeros(size(y));

        oscillatory_forcing = w*A*cos(w*t);

        V = V_0 + V_1.*cos(w*t - delta);
        
        theta = y(1);
        psi = y(2);
        phi = y(3);
        
        f1 = -B*(sin(2*theta)*sin(2*phi))/4;
        f2 = (B/2)*cos(theta)*cos(2*phi);
        f3 = (1 - B*cos(2*phi))/2;
        
        % Rotational dynamics
        dy(1) = oscillatory_forcing*cos(psi) + G*f1;
        dy(2) = - oscillatory_forcing.*cot(theta).*sin(psi) + G*f2;
        dy(3) = oscillatory_forcing.*sin(psi)./sin(theta) + G*f3;

        % Translational dynamics
        dy(4) = V(1)*cos(theta) + V(2)*sin(psi)*sin(theta) + V(3)*cos(psi)*sin(theta);
        dy(5) = V(1)*sin(phi)*sin(theta) + V(2)*(cos(phi)*cos(psi) - sin(phi)*cos(theta)*sin(psi)) ...
            - V(3)*(cos(phi)*sin(psi) + cos(theta)*sin(phi)*cos(psi));
        dy(6) = -V(1)*(cos(phi)*sin(theta)) + V(2)*(sin(phi)*cos(psi) + cos(theta)*cos(phi)*sin(psi))...
            + V(3)*(cos(theta)*cos(phi)*cos(psi) - sin(phi)*sin(psi)) + G*y(5);
        
    end

end

