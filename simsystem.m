function [yhat,xhat_out] = simsystem(A,B,C,D,x0,u)
%% Instructions:
% Implement a function that simulates the system here!
% Use the following function inputs and outputs.

% Function INPUT 
% A         System matrix A (matrix of size n x n)
% B         System matrix B (matrix of size n x m)
% C         System matrix C (matrix of size l x n)
% D         System matrix D (matrix of size l x m)
% x0        Initial state (vector of size n x one)
% u         system input (matrix of size N x m)

% Function OUTPUT
% yhat      predicted output (vector of size l x one)

% Initialize matrices
sim_len = size(u,1);
xhat = zeros(size(A,1),sim_len+1);
yhat = zeros(sim_len,1);

% Setting initial point for xhat
xhat(:,1) = x0;

% Simulation
for k=1:sim_len
    xhat(:,k+1) = A*xhat(:,k) + B*u(k,:)';
    yhat(k) = C*xhat(:,k) + D*u(k,:)';
end

% Setting second output
xhat_out = xhat(:,1:sim_len)';

end