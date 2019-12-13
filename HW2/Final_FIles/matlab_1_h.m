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
A = 5000; omega_x=5; omega_y=7;

% Define phi_h and phi_L. Because -1=<cos<=1 in the specified region, 
%phi(x,y) can get values between -A and A
phi_L = -A;
phi_H = A;

%calculate the function values 
phi_xy=A*cos(2*pi*omega_x*(x_grid)).* cos(2*pi* omega_y *(y_grid));

%Plot the function as an image histogram.
imshow(phi_xy, [phi_L phi_H]);
title('phi(x,y) image');

figure('Name', 'phi(x,y)');
surf(x_grid,y_grid,phi_xy, 'LineStyle' , 'none');

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
% energy_der_phi_y = 1.12113 * 10^10

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

% Parameters we got: 
%Nx = 36.6292
%Ny = 51.2788
%b = 2.662
%MSE_opt = 9.7577 * 10^5
% Do same optimization for B_high
con_h = @(x)deal([], (x(1)*x(2)*x(3)-B_high));
lb_h = [0.001 0.001 0.001];
ub_h = [B_high B_high B_high];
[opt_h, opt_mse_h] = fmincon(MSE, x0, [], [], [],[], lb_h, ub_h, con_h, options);
% Parameters we got:
%Nx = 90.8949
%Ny = 127.2479
%b = 4.3229
%MSE_opt = 1.4548*10^5

%% Part f + g

% Hold the parameters result
best_parameters_low = zeros(1,3);

% Hold best MSE result
best_mse_low = intmax;

% Run over all possible values for Nx
for Nx = 1 : B_low
    
    % Run over all possibles values for Ny
   for Ny = 1 : B_low
       
      % Calculate b accordiny to current Nx and Ny
      % If b=0, this run in useless, continue
      b = floor(B_low/(Nx*Ny));
      if b==0
          continue;
      end
      x = [Nx, Ny, b];
      
      % Calculate the MSE, if its better than the best one until this run,
      % update the parameters and save this result
      curr_mse = MSE(x);
      if(curr_mse < best_mse_low)
          best_mse_low = curr_mse;
          best_parameters_low = [Nx Ny b];
      end    
   end  
end

% Parameters we got foo B_low:
% Nx = 34
% Ny = 49
% b =3
% MSE_low = 9.9613*10^5

% Do the same for B_high...

best_parameters_high = zeros(1,3);
best_mse_high = intmax;
for Nx = 1 : B_high
   for Ny = 1 : B_high
      b = floor(B_high/(Nx*Ny));
      if b==0
          continue;
      end
      x = [Nx, Ny, b];
      curr_mse = MSE(x);
      if(curr_mse < best_mse_high)
          best_mse_high = curr_mse;
          best_parameters_high = [Nx Ny b];
      end    
   end  
end

% Parameters we got for B_high:
% Nx = 96
% Ny = 130
% b =4
% MSE_low = 7.0049*10^4


% Plot the results: 

% First, plot for values obtained by B_low

continuous_approx_delta_x = 1/best_parameters_low(1);
continuous_approx_delta_y = 1/best_parameters_low(2);
%2D grid of (0,1]x(0,1]
grid_x = 0 : continuous_approx_delta_x : 1;
grid_y = 0 : continuous_approx_delta_y : 1;
[x_grid,y_grid] = meshgrid(grid_x,grid_y);

%calculate the function values 
phi_xy_l=A*cos(2*pi*omega_x*(x_grid)).* cos(2*pi* omega_y *(y_grid));

%Plot the function as an image histogram.
imshow(phi_xy_l, [phi_L phi_H]);
title('phi(x,y) image reconstructed with B low');

figure('Name', 'phi(x,y) reconstructed with B low');
surf(x_grid,y_grid,phi_xy_l, 'LineStyle' , 'none');
title('phi(x,y) reconstructed with B low');

% Do the same for B_high

continuous_approx_delta_x = 1/best_parameters_high(1);
continuous_approx_delta_y = 1/best_parameters_high(2);
%2D grid of (0,1]x(0,1]
grid_x = 0 : continuous_approx_delta_x : 1;
grid_y = 0 : continuous_approx_delta_y : 1;
[x_grid,y_grid] = meshgrid(grid_x,grid_y);

%calculate the function values 
phi_xy_h=A*cos(2*pi*omega_x*(x_grid)).* cos(2*pi* omega_y *(y_grid));

%Plot the function as an image histogram.
imshow(phi_xy_h, [phi_L phi_H]);
title('phi(x,y) image reconstructed with B high');

figure('Name', 'phi(x,y) reconstructed with B high');
surf(x_grid,y_grid,phi_xy_h, 'LineStyle' , 'none');
title('phi(x,y) reconstructed with B high');

