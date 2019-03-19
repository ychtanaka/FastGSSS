function [selected_nodes, T] = fastGSSS(L_n, F, bw, nu, eigen_decomp, varargin)
% FASTGSSS Eigendecomposition-free sampling set selection for graphs
%
% Usage:  [selected_nodes, T_g_tmp1] = fastGSSS(L_n, F, eigen_decomp, U, E, m, bw, alpha)
% Example: See demo
%
% Input parameter (G is required, others are optional):
% L_n: Symmetric normalized graph Laplacian.
% F: # of vertices to be selected.
% bw: (Estimated) bandwidth the graph signals.
% nu: parameter controlling the width of the filter kernel.
% eigen_decomp: flag for eigendecomposition.
%     1: Performs exact filtering (needs U and E)
%     0: Approximates filter with Chebyshev polynomial approximation (needs m)
% varargin: Inputs required for sensor selection.
%     If eigen_decomp = 1, varargin{1} = U and varargin{2} = E
%       U: Eigenvector matrix of L_n. i.e., GFT matrix.
%       E: Eigenvalue matrix of L_n
%     If eigen_decomp = 0, varargin{1} = m
%       m: Polynomial order.
% 
% Output parameters:
% selected_nodes: Indices of vertex selected.
% T: Basis for reconstruction.
% 
% References:
% A. Sakiyama, Y. Tanaka, T. Tanaka, and A. Ortega,
% "Eigendecomposition-free sampling set selection for graph signals," IEEE
% Transactions on Signal Processing, accepted.
%
% Copyright (c) 2019 Akie Sakiyama, Yuichi Tanaka, Toshihisa Tanaka, and
% Antonio Ortega

%% Input check
if eigen_decomp == 1

    if nargin ~= 7
        error('For eigen_decomp = 1, U and E must be declared!');
    else
        U = varargin{1};
        E = varargin{2};
    end

else
    
    if nargin ~= 6
        error('For eigen_decomp = 0, m must be declared!');
    else
        m = varargin{1};
    end
    
end

%% parameters
N=length(L_n);
A_tmp = logical(L_n - diag(diag(L_n)));
numedge = sum(A_tmp(:))/2;
p = numedge/N; % edge probability
n = F/N; % sampling ratio
k = bw/N; %bandwidth

%% Preparing T
if eigen_decomp==1
    
    lmax=max(E);
    g_E=exp(-nu*p*n*k*E/lmax);
    T_g_tmp1=U*diag(g_E)*U.';

else
    arange=[0,2];
    g = @(x)(exp(-nu*p*n*k*x/2));
    c = sgwt_cheby_coeff(g, m, m+1, arange);
    T_g_tmp1 = sgwt_cheby_op2(L_n, c, arange);
end

%% SSS
T_g_tmp=abs(T_g_tmp1);
sensor = selection([],T_g_tmp);
selected_nodes = sensor;

for i=1:F-1
    sensor = selection(selected_nodes,T_g_tmp);
    selected_nodes = [selected_nodes sensor];
end

T = T_g_tmp1;

end

%% Internal SSS function
function sensor = selection(selected,T_g_tmp)

if ~isempty(selected)
    T=sum(T_g_tmp(:,selected),2);
    T2=mean(T)-T;
    T2(logical(T2<0))=0;
    T_g=(T_g_tmp)*(T2);   
else
     T_g=sum(T_g_tmp);
end
    
T_g(selected) = 0;
[~,sensor]  =max(T_g);

end

