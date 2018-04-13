%% -- load all files
addpath(genpath('./bin'));
addpath(genpath('./class'));
addpath(genpath('./contrib'));
addpath(genpath('./femm'));


%% -- create Aniso class
an = aniso();
nodes = an.rte.getNodes();

% nodes are aligned in column major.
% the diffusion preconditioner can be created on the coarser grid.
%% -- set the integral
ss = @(x)(50);
aa = @(x)(0.2);

an.diffgen(ss, aa);
an.setCoeff(ss, aa);

%% set the preconditioner.

%% gmres
chargeFun = @(x) (exp(- 25 * ((x(:,1)-0.5).^2 + (x(:,2)-0.5).^2)));
charge = chargeFun(nodes);


%%
[y, flag, relres, iter, resvec] = an.solve(charge);

