addpath(genpath('./bin'));
addpath(genpath('./class'));
addpath(genpath('./contrib'));

rte = Aniso(64, 3, 1, 0.8, 8, 4, 20);
nodes = rte.getNodes();
n = size(nodes, 1);

sigma_s = 10*ones(n,1);
sigma_t = 0.2 + sigma_s;

rte.setCoeff(sigma_s, sigma_t);

tic;
rte.cache();
toc;
%%
chargeFun = @(x) (exp(- 5 * ((x(:,1)-0.5).^2 + (x(:,2)-0.5).^2)));

charge = chargeFun(nodes);
rhs = rte.mapping(charge);
%%
A = @(x)(x - sigma_s .* rte.mapping(x));
tic;
[y, flag, relres, iter, resvec] = gmres(A, rhs, 40, 1e-12, 400);
toc;

