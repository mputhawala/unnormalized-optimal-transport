mex -output MEXUnbalancedOTSolver1D CXXFLAGS="\$CXXFLAGS -std=c++11" MexUnnormalizedOtSolver1DEntry.cpp UnnormalizedOTSolver1D.cpp FunctionUtilities.cpp% -D_MATLAB_DEBUG

%% Test the build
max_iters = 2000;
n_x = 35;
n_t = 10;
tau_1 = 1e-3;
tau_2 = 1e-1;
alpha = 1;
enforce_zero_f = false;
obj_stopping_val = -1;
constraint_stopping_val = -1;

%% gaussian params
mu_0 = 0;
sigma_0 = 1e-1;
mu_1 = 2/3;
sigma_1 = 1e-1;


gaussian = @(x, mu, sigma) exp(-(x - mu).^2 ./ (2 * sigma.^2));
ind2x = @(ind) (ind - 1) / (n_x - 1);

rho_0 = zeros([n_x, 1]);
rho_1 = zeros([n_x, 1]);

for i = 1:n_x
    x = ind2x(i);
    %[x,y]
    rho_0(i) = gaussian(x, mu_0, sigma_0);
    rho_1(i) = gaussian(x, mu_1, sigma_1);
end

x_cell_width = 1/(n_x - 1);
rho_0 = rho_0/(sum(rho_0(:))*x_cell_width);
rho_1 = rho_1/(sum(rho_1(:))*x_cell_width);

%%

[emd, m_x, phi, u, f, misc] = MEXUnbalancedOTSolver(rho_0, rho_1, n_t, tau_1, tau_2, alpha, enforce_zero_f, obj_stopping_val, constraint_stopping_val, max_iters);