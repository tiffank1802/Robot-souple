function [B] = Bconstruct(M,C,K,nA,dt,beta,gamma)
nx=size(M,1);
B = zeros(3*nx, nx);
for j = 1:nx
    u0=zeros(nx,1);
    v0=zeros(nx,1);
    a0=zeros(nx,1);
    eta=zeros(nx,1);
    eta(j)=1;
    [u,v,a]=newmark1stepMRHS(M,C,K,eta,u0,v0,a0,dt,beta,gamma);
    B(:,j) = [u;v;a];
end
end