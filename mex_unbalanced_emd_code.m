mex -output MEXUnbalancedOTSolver2D CXXFLAGS="\$CXXFLAGS -std=c++11" main.cpp UnnormalizedOTSolver.cpp FunctionUtilities.cpp

%% Test the build
max_iters = 2000;
n_x = 15;
n_y = 15;
n_t = 5;
tau_1 = 2e-3;
tau_2 = 1e-1;
alpha = 1;
enforce_zero_f = false;
obj_stopping_val = -1;
constraint_stopping_val = -1;

%% gaussian params
mu_0 = [0, 0];
sigma_0 = 1e-1*[1, 1];
mu_1 = [1, 1];
sigma_1 = 1e-1*[1, 1];

%ar(i, j) = exp(-pow(x - mu_1, 2) / (2.0 * pow(sigma_1, 2)) - pow(y - mu_2, 2) / (2.0 * pow(sigma_2, 2)));

gaussian = @(x, mu, sigma) prod(exp(-(x - mu).^2 ./ (2 * sigma.^2)));
ind2x = @(ind) (ind - 1) / (n_x - 1);
ind2y = @(ind) (ind - 1) / (n_y - 1);

rho_0 = zeros([n_x, n_y]);
rho_1 = zeros([n_x, n_y]);

for i = 1:n_x
    for j = 1:n_y
        x = ind2x(i);
        y = ind2y(j);
        %[x,y]
        rho_0(i, j) = gaussian([x, y], mu_0, sigma_0);
        rho_1(i, j) = gaussian([x, y], mu_1, sigma_1);
    end
end

x_cell_width = 1/(n_x - 1);
y_cell_width = 1/(n_y - 1);
rho_0 = rho_0/(sum(rho_0(:))*x_cell_width*y_cell_width);
rho_1 = rho_1/(sum(rho_1(:))*x_cell_width*y_cell_width);

%%

[emd, m_x, m_y, phi, u, f, misc] = MEXUnbalancedOTSolver(rho_0, rho_1, n_t, tau_1, tau_2, alpha, enforce_zero_f, obj_stopping_val, constraint_stopping_val, max_iters);