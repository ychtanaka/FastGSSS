function mem_fn_recon = reconstruction_fastGSSS(mem_fn,S_opt,T,k)
%RECONSTRUCTION_POCS Summary of this function goes here
%   Detailed explanation goes here

S_set=find(S_opt);
T_k=T^k;
mem_fn_recon=T_k(:,S_set)*((T_k(S_set,S_set))\mem_fn(S_set,:));


