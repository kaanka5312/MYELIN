function [x, time] = wilsoncowan_RK2(tau, b, W, k, s, C, tspan, dt, burn)
arguments
    tau                    % timescale
    b                      % bias term
    W                      % connectivity matrix
    k                      % sharpness of f-I curve
    s                      % input
    C                      % coupling strength
    tspan = 600            % time (s)
    dt = 10 / 1000         % time step (s)
    burn = 100             % burnin period (s)                
end

nrois = length(W);
time = 0:dt:tspan;
ntime = length(time);
s_input = s + randn(nrois, ntime);
x = zeros(nrois, ntime);
x(:, 1) = randn(nrois, 1);

for t = 1:(ntime-1)
    k1 =(dt / 2) * wilsoncowan_dxdt(x(:, t), ...
        s_input(:, t), tau, b, W, k, C);
    k2 = dt * wilsoncowan_dxdt(x(:, t) + k1, ...
        s_input(:, t), tau, b, W, k, C);
    x(:, t + 1) = x(:, t) + k2;
end

burn_samples = round(burn / dt);
x = x(:, burn_samples:end);
time = time(:, burn_samples:end);
end

function y = f(k, x)
y = 1 ./ (1 + exp(-k * x));
end

function dxdt = wilsoncowan_dxdt(x, s_input, tau, b, W, k, C)
dxdt = (-x + f(k, C .* W * x + b + s_input)) ./ tau;
end