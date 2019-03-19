function [selected_nodes, T, G] = gsp_demo_fastGSSS(gname, numnodes)
% GSP_DEMO_FASTGSSS Demonstration of eigendecomposition-free sampling set
% selection for graph signals
%
% This demonstration performs eigendecomposition-free sampling set
% selection for graph signals.
%
% Usage:  [selected_nodes, T, G] = gsp_demo_fastGSSS(gname, numnodes);
%
% Input parameter (optional):
% gname: Name of the graph. For visualization, the graph must have its
% vertex coordinates.
% numnodes: # of vertices to be selected.
% 
%
% Output parameters:
% selected_nodes: Indices for selected vertices.
% T: Basis for reconstruction.
% G: Graph.
%
% References:
% A. Sakiyama, Y. Tanaka, T. Tanaka, and A. Ortega,
% "Eigendecomposition-free sampling set selection for graph signals," IEEE
% Transactions on Signal Processing, accepted.
%
% Copyright (c) 2019 Akie Sakiyama, Yuichi Tanaka, Toshihisa Tanaka, and
% Antonio Ortega

%% input check
if nargin == 0
    gname = 'community';
    numnodes = 10;
elseif nargin == 1
    numnodes = 10;
else
end

switch gname
    case 'community'
        G = gsp_community(500);
    case 'sensor'
        param_graph.connected = 1;
        param_graph.distribute = 1;
        param_graph.N_try = 100;
        G = gsp_sensor(500, param_graph);
end

%% preparing parameters and graph (signal)
order = 12; % polynomial approximation order
nu = 75; % parameter controlling the width of the filter kernel
bw = 100; % bandwidth
N = G.N;

A = G.W;
L = diag(sum(A,2))-A;
L_n = diag(diag(L).^(-0.5))*L*diag(diag(L).^(-0.5));

[eigenvectors,eigenvalues,~] = svd(full(L_n));

% Sort eigenvectors and eigenvalues
[E, inds] = sort(diag(eigenvalues),'ascend');
eigenvectors=eigenvectors(:,inds');
% Set first component of each eigenvector to be nonnegative
signs=sign(eigenvectors(1,:));
signs(signs==0)=1;
U = eigenvectors*diag(signs);

% graph signal
f = U*([sqrt(0.2)*randn(bw,1);zeros(G.N-bw,1)]+ sqrt(5* 10^(-3))*randn(length(G.N),1));
G.L = L_n;
f_or = f;

%% SSS
tic;
[selected_nodes, T] = fastGSSS(L_n, numnodes, bw, nu, 0, order);
tElasped = toc;
fprintf('Execution time for selecting %d (out of %d) vertices: %f seconds \n', numnodes, N, tElasped);

sampling_set = zeros(N,1);
sampling_set(selected_nodes) = 1;

% reconstruction
f_out = reconstruction_fastGSSS(f_or, sampling_set, T, order);

MSE = mean((f_or - f_out).^2);
fprintf('%d vertices are selected from %d vertices: MSE = %g \n', numnodes, N, MSE);

% visualization
samp_vertices = -1 * ones(G.N, 1);
samp_vertices(selected_nodes) = 1;

figure('Color',[1 1 1]);
param.vertex_size = 200;
param.show_edges = 1;
gsp_plot_signal(G, samp_vertices, param);
title('Sampled vertices');

