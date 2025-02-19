function time_plots
%TIME_PLOTS Plot four graphs on one figure, showing sin(theta)*sin(phi),
%sin(theta)*cos(phi), sin(theta)*sin(psi), and sin(theta)*cos(psi), and
%another three graphs for the position on another figure.

% G = gamma: the shear rate of the flow
G = 1;

% A: the magnitude of yawing
A = 2;

% B: the Bretheton coefficient
B = 0.9;

% w: the yawing frequency.
w = 3;

% V_i is the velocity in the \hat{e}_i direction (i.e. in a direction fixed
% in the swimmer frame). It has the form V_i = V_0i + V_1i*cos(w*t -
% delta_i).

V_01 = -0.2;
V_02 = 0.5;
V_03 = 0.2;

V_11 = 0.2;
V_12 = 0.6;
V_13 = 0.5;

delta_1 = pi/2;
delta_2 = pi/4;
delta_3 = -pi/4;


V_0 = [V_01,V_02,V_03];
V_1 = [V_11,V_12,V_13];
delta = [delta_1, delta_2, delta_3];


% Initial conditions for the orientation via theta, phi, and psi.
Init_theta = pi/6;
Init_psi = pi/12;
Init_phi = pi/12;

% The initial conditions vector for full simulations. We set the initial
% position to be the origin without loss of generality.
Init = [Init_theta,Init_psi,Init_phi, 0, 0, 0];

% Time interval
T_init = 0;
T_max = 100;

% ODE solver tolerance
options = odeset('RelTol',3e-14,'AbsTol',3e-14);

% Solve full ODEs
tic
[t,y] = solve_ODE(G,A,B,V_0,V_1,delta,w,T_init,T_max,Init,options);
toc

% Relabel numerical solutions
Theta = y(:,1);
Psi = y(:,2);
Phi = y(:,3);
X = y(:,4);
Y = y(:,5);
Z = y(:,6);


% Solve full ODEs with zero flow
tic
[t_no_flow,y_no_flow] = solve_ODE(0,A,B,V_0,V_1,delta,w,T_init,T_max,Init,options);
toc

% Relabel numerical solutions
Theta_no_flow = y_no_flow(:,1);
Psi_no_flow = y_no_flow(:,2);
Phi_no_flow = y_no_flow(:,3);
X_no_flow = y_no_flow(:,4);
Y_no_flow = y_no_flow(:,5);
Z_no_flow = y_no_flow(:,6);


% Precalculate Bessel functions for slow evolution equations
J0_1A = besselj(0,A);
J1_1A = besselj(1,A);
J0_2A = besselj(0,2*A);


% Solve slow evolution ODEs
t_span = t;
tic
[t_slow,y_slow] = solve_slow_ODE(G,B,V_0,V_1,delta,J0_1A,J1_1A,J0_2A,t_span,Init,options);
toc


% Relabel slow evolution solutions
alpha = y_slow(:,1);
beta = y_slow(:,2);
gamma = y_slow(:,3);
X_asy = y_slow(:,4);
Y_asy = y_slow(:,5);
Z_asy = y_slow(:,6);

% Define line width and font size for figures
L = 4;
F = 26;

% Specify yawing function
f = A*sin(w*t);

% Define the asymptotic result for cos(theta)
cos_theta_asy = cos(alpha).*cos(f) - sin(alpha).*cos(beta).*sin(f);

% Colour for figures
blue = [0.30,0.75,0.93];

% Set up figure
figure('units','normalized','outerposition',[0 0 1 1],'Renderer','painters')
Tile = tiledlayout(4,1);

% Define the asymptotic result for sin(theta)*cos(phi)
sin_theta_cos_phi_asy = (cos(alpha).*cos(beta).*sin(f) + sin(alpha).*cos(f)).*cos(gamma) ...
    - sin(beta).*sin(f).*sin(gamma);

% Plot the results for sin(theta)*cos(phi)
ax1 = nexttile;
plot(t,sin(Theta).*cos(Phi),'Color',blue,'LineWidth',L+1)
hold on
plot(t_no_flow,sin(Theta_no_flow).*cos(Phi_no_flow),'r','LineWidth',L-2)
plot(t_slow,sin_theta_cos_phi_asy,':','Color','k','LineWidth',L)

ylim([-1 1])
% xlabel('$t$','Interpreter','latex','FontSize',F);
ylabel('$\sin \theta \cos \phi$','Interpreter','latex','FontSize',F);
set(gcf,'Color',[1,1,1])
set(gca,'TickLabelInterpreter','latex','FontSize',F);

% Define the asymptotic result for sin(theta)*sin(phi)
sin_theta_sin_phi_asy = (cos(alpha).*cos(beta).*sin(f) + sin(alpha).*cos(f)).*sin(gamma) ...
    + sin(beta).*sin(f).*cos(gamma);

% Plot the results for sin(theta)*sin(phi)
ax2 = nexttile;
plot(t,sin(Theta).*sin(Phi),'Color',blue,'LineWidth',L+1)
hold on
plot(t_no_flow,sin(Theta_no_flow).*sin(Phi_no_flow),'r','LineWidth',L-2)
plot(t_slow,sin_theta_sin_phi_asy,':','Color','k','LineWidth',L)

ylim([-1 1])
% xlabel('$t$','Interpreter','latex','FontSize',F);
ylabel('$\sin \theta \sin \phi$','Interpreter','latex','FontSize',F);
set(gcf,'Color',[1,1,1])
set(gca,'TickLabelInterpreter','latex','FontSize',F);

% Define the asymptotic result for sin(theta)*cos(psi) and
% sin(theta)*sin(psi)
sin_theta_cos_psi_asy = cos(alpha).*sin(f) + sin(alpha).*cos(beta).*cos(f);
sin_theta_sin_psi_asy = sin(alpha).*sin(beta);

% Plot the results for sin(theta)*cos(psi)
ax3 = nexttile;
plot(t,sin(Theta).*cos(Psi),'Color',blue,'LineWidth',L+1)
hold on
plot(t_no_flow,sin(Theta_no_flow).*cos(Psi_no_flow),'r','LineWidth',L-2)
plot(t,sin_theta_cos_psi_asy,':','Color','k','LineWidth',L)

ylim([-1 1])
% xlabel('$t$','Interpreter','latex','FontSize',F);
ylabel('$\sin \theta \cos \psi$','Interpreter','latex','FontSize',F);
set(gcf,'Color',[1,1,1])
set(gca,'TickLabelInterpreter','latex','FontSize',F);

% Plot the results for sin(theta)*sin(psi)
ax4 = nexttile;
plot(t,sin(Theta).*sin(Psi),'Color',blue,'LineWidth',L+1)
hold on
plot(t_no_flow,sin(Theta_no_flow).*sin(Psi_no_flow),'r','LineWidth',L-2)
plot(t,sin_theta_sin_psi_asy,':','Color','k','LineWidth',L)

ylim([-1 1])
% yticks([-1 0])
% xlabel('$t$','Interpreter','latex','FontSize',F);
ylabel('$\sin \theta \sin \psi$','Interpreter','latex','FontSize',F);
set(gcf,'Color',[1,1,1])
set(gca,'TickLabelInterpreter','latex','FontSize',F);

legend('Numerical','Ignoring slow evolution','Asymptotic','Interpreter','latex','FontSize',F,...
    'Orientation','Horizontal','Location','SouthEast');



linkaxes([ax1,ax2,ax3,ax4],'x');
xlabel(Tile,'$t$','Interpreter','latex','FontSize',F);

xticklabels(ax1,{})
xticklabels(ax2,{})
xticklabels(ax3,{})
Tile.Padding = 'tight';
Tile.TileSpacing = 'tight';

set(gcf,'Color',[1,1,1])
set(gca,'TickLabelInterpreter','latex','FontSize',F);

%%% Save figure
% FIG = gcf;
% exportgraphics(FIG,'Orientation_over_time.eps')




% Generate figure for translational dynamics
figure('units','normalized','outerposition',[0 0 1 1],'Renderer','painters')
Tile = tiledlayout(3,1);

% Plot the results for X
ax1 = nexttile;
plot(t,X,'Color',blue,'LineWidth',L+1)
hold on
plot(t_no_flow,X_no_flow,'r','LineWidth',L-2)
plot(t_slow,X_asy,':','Color','k','LineWidth',L)
% xlabel('$t$','Interpreter','latex','FontSize',F);
ylabel('$X$','Interpreter','latex','FontSize',F);
set(gcf,'Color',[1,1,1])
set(gca,'TickLabelInterpreter','latex','FontSize',F);

% Plot the results for Y
ax2 = nexttile;
plot(t,Y,'Color',blue,'LineWidth',L+1)
hold on
plot(t_no_flow,Y_no_flow,'r','LineWidth',L-2)
plot(t_slow,Y_asy,':','Color','k','LineWidth',L)
% xlabel('$t$','Interpreter','latex','FontSize',F);
ylabel('$Y$','Interpreter','latex','FontSize',F);
set(gcf,'Color',[1,1,1])
set(gca,'TickLabelInterpreter','latex','FontSize',F);

% Plot the results for Z
ax2 = nexttile;
plot(t,Z,'Color',blue,'LineWidth',L+1)
hold on
plot(t_no_flow,Z_no_flow,'r','LineWidth',L-2)
plot(t_slow,Z_asy,':','Color','k','LineWidth',L)

% xlabel('$t$','Interpreter','latex','FontSize',F);
ylabel('$Z$','Interpreter','latex','FontSize',F);
set(gcf,'Color',[1,1,1])
set(gca,'TickLabelInterpreter','latex','FontSize',F);



legend('Numerical','Ignoring slow evolution','Asymptotic','Interpreter','latex','FontSize',F,...
    'Orientation','Horizontal','Location','SouthWest');



linkaxes([ax1,ax2,ax3],'x');
xlabel(Tile,'$t$','Interpreter','latex','FontSize',F);

xticklabels(ax1,{})
xticklabels(ax2,{})
Tile.Padding = 'tight';
Tile.TileSpacing = 'tight';

set(gcf,'Color',[1,1,1])
set(gca,'TickLabelInterpreter','latex','FontSize',F);

%%% Save figure
% FIG = gcf;
% exportgraphics(FIG,'Position_over_time.eps')


end

