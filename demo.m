%% -- load all files
addpath(genpath('./bin'));
addpath(genpath('./class'));
addpath(genpath('./contrib'));
addpath(genpath('./femm'));


%% -- create Aniso class
an = aniso(0.8, 5);
nodes = an.rte.getNodes();

% nodes are aligned in column major.
% the diffusion preconditioner can be created on the coarser grid.
%% -- set the integral
ss = @(x)(20);
aa = @(x)(0.2);

an.diffgen(ss, aa);
an.setCoeff(ss, aa);

%% set the preconditioner.

%% gmres
chargeFun = @(x) (exp(- 25 * ((x(:,1)-0.5).^2 + (x(:,2)-0.5).^2)));
charge = chargeFun(nodes);

n = size(charge, 1);
Charge = zeros(an.N * n, 1);
Charge(1:an.n) = charge;

%%
[y, flag, relres, iter, resvec] = an.solve(Charge, 0);

