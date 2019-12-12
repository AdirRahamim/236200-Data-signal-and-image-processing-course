%% Matlab - part I code

clear; clc;
%% section b

% We want to approximate the signal in a very high resolution,
% we choose delta between each sample to be very small.

continuous_approx_delta = 0.001;
%2D grid of (0,1]x(0,1]
grid = 0 : continuous_approx_delta : 1;
[x_grid,y_grid] = meshgrid(grid,grid);

%function parameters
A = 5000; omega_x=5; omega_y=3;

% Define phi_h and phi_L. Because -1=<cos<=1 in the specified region, 
%phi(x,y) can get values between -A and A
phi_L = -A;
phi_H = A;

%calculate the function values 
phi_xy=A*cos(2*pi*omega_x*(x_grid)).* cos(2*pi* omega_y *(y_grid));

%Plot the function as an image histogram.
%imshow(phi_xy, [phi_L phi_H]);
%title('phi(x,y) image');

%figure('Name', 'phi(x,y)');
%surf(x_grid,y_grid,phi_xy, 'LineStyle' , 'none');

%% section c

% Calulate numerically phi_L and phi_H by finding the minimum and maximum
% element of the phi_xy matrix

phi_L = min(min(phi_xy));
phi_H = max(max(phi_xy));
value_range = phi_H - phi_L;

% Calculate the derivative of phi according to x the derivative
% according to y in every point in the area. We will calulate it using
% the derivative definistion. 
% First we calculate (x+delta_x) and (y+delta_y)
phi_xy_delta_x = A*cos(2*pi*omega_x*(x_grid+continuous_approx_delta)).* cos(2*pi* omega_y *(y_grid));
phi_xy_delta_y = A*cos(2*pi*omega_x*(x_grid)).* cos(2*pi* omega_y *(y_grid+continuous_approx_delta));

% Now we calculate approximation for the derivative in every point by
% Fx = (F(x+delta_x)-F(x))/delta_x
der_phi_x = (phi_xy_delta_x-phi_xy)./continuous_approx_delta;
der_phi_y = (phi_xy_delta_y-phi_xy)./continuous_approx_delta;

% Lastly, we will aproximate the integral of the energy of the 
% parial drivatives of phi, we will approximate it by calculate
% the sum of the areas of the 3D squres created beneath the 
% energy functions

energy_der_phi_x = sum(sum((der_phi_x.^2).*((continuous_approx_delta).^2)));
energy_der_phi_y = sum(sum((der_phi_y.^2).*((continuous_approx_delta).^2)));

% The values we got: 
% Values range = 2A = 10,000
% energy_der_phi_x = 6.1803 * 10^9
% energy_der_phi_y = 2.225 * 10^9

%% Section d + e

% Using the numerical estimations from section c we need to determine 
% Nx, Ny, b that satisfy the nonlinear optimization problem:
% min (MSE_TOT(Nx, Ny, b)) s.t Nx * Ny * b = B

% Define the possible bit-budgets
B_low = 5000;
B_high = 50000;

% Define the MSE function as we seen in class.
% x(1) represents Nx, x(2) represents Ny and x(3) represents b
MSE = @(x)((1/(12*(x(1))^2))*energy_der_phi_x) + ...
    ((1/(12*(x(2))^2))*energy_der_phi_y) + ...
    (((value_range)^2)/(12*(2^(2*x(3)))));

% Define starting values for x for finding the best optimization.
% We dont set them zero that fmincon wont fail because 
% division by 0
x0 = [0.001 0.001 0.001]';

% Define constraits. We have only equality constraint
% Nx*Ny*b-B=0.
con_l = @(x)deal([], (x(1)*x(2)*x(3)-B_low));

% Define lower and upper bounds for the values of 
% Nx, Ny, b
lb_l = [0.001 0.001 0.001];
ub_l = [B_low B_low B_low];

% We need to set options to increase max number of iterations
options = optimoptions('fmincon', 'MaxFunctionEvaluation', inf, 'MaxIterations' , inf);
[opt_l, opt_mse_l] = fmincon(MSE, x0, [], [], [],[], lb_l, ub_l, con_l, options);

% Do same optimization for B_high
con_h = @(x)deal([], (x(1)*x(2)*x(3)-B_high));
lb_h = [0.001 0.001 0.001];
ub_h = [B_high B_high B_high];
[opt_h, opt_mse_h] = fmincon(MSE, x0, [], [], [],[], lb_h, ub_h, con_h, options);


%% Part f + g

%Set bit budget
B = B_low;

best_parameters = zeros(1,3);
best_mse = intmax;
tic;
for Nx = 1 : B
   for Ny = 1 : B
      b = floor(B/(Nx*Ny));
      if b==0
          continue;
      end
      %fprintf('Been here\n');
      x = [Nx, Ny, b];
      curr_mse = MSE(x);
      if(curr_mse < best_mse)
          best_mse = curr_mse;
          best_parameters = [Nx Ny b];
      end
      
   end
    
end
toc;
