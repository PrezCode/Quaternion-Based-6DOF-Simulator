%Quaternion State Plot
clf;
states = ['Object_Simulation_Data.txt'];
data = readmatrix(states);

time = data(:, 1);
X = data(:, 2);
Y = data(:, 3);
Z = data(:, 4);
U = data(:, 5);
V = data(:, 6);
W = data(:, 7);
P = data(:, 8);
Q = data(:, 9);
R = data(:, 10);
phi = data(:, 11);
theta = data(:, 12);
psi = data(:, 13);

% Displacement Plot
figure(1);
plot(time, X/1000, '--','Color',"r"); 
hold on; 
plot(time, Y/1000, '--','Color',"g");
plot(time, -Z/1000,'--','Color',"b");
hold off;
grid on;
xlabel('Seconds');
ylabel('Kilometers');
title('Displacement over Time');
legend('X','Y','Z');
% Velocity Plot
figure(2);
plot(time, U, '--','Color',"r"); 
hold on; 
plot(time, V, '--','Color',"g");
plot(time, W,'--','Color',"b");
hold off;
grid on;
xlabel('Seconds');
ylabel('Meters per Second');
title('Velocity over Time');
legend('U_e','V_e','W_e');
% Rotation Rate Plot
figure(3);
plot(time, P*180.0/pi, '--','Color',"r"); 
hold on; 
plot(time, Q*180.0/pi, '--','Color',"g");
plot(time, R*180.0/pi,'--','Color',"b");
hold off;
grid on;
xlabel('Seconds');
ylabel('Degrees per Second');
title('Angular Velocity over Time');
legend('Roll(P)','Pitch(Q)','Yaw(R)');
% Euler Angle Plot
figure(4);
plot(time, phi*180.0/pi, '--','Color',"r"); 
hold on; 
plot(time, theta*180.0/pi, '--','Color',"g");
plot(time, psi*180.0/pi,'--','Color',"b");
hold off;
grid on;
xlabel('Seconds');
ylabel('Degrees');
title('Euler Angles over Time');
legend('Roll(phi)','Pitch(theta)','Yaw(psi)');
%3-D Displacement Plot
figure(5);
plot3(X/1000, Y/1000, Z/1000, 'LineWidth', 3);
hold on;
radius = 6371.0082;
[xsphere, ysphere, zsphere] = sphere(100);
surf(xsphere*radius, ysphere*radius, zsphere*radius);
plot3(0, 0, 0);
hold off;
grid on;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('Position over Time');