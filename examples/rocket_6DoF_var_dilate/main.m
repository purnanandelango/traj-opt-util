clearvars
clc

prb = problem_data(30,10,2e2,0.1,0.1);

load('recent_solution','x','u','tau');
[xbar,ubar] = misc.create_initialization(prb,1,x,u,tau);

[xbar,ubar] = scp.run_ptr_noparam_varscl(xbar,ubar,prb,@sys_cnstr_cost);
tvecbar = prb.time_grid(prb.tau,xbar,ubar);

% Simulate solution on fine grid

% Simulate on [0,1] grid
[tau,x,u] = disc.simulate_dyn(xbar(:,1),{prb.tau,ubar},@(t,x,u) prb.dyn_func(t,x,u),[0,1],prb.Kfine,prb.disc);
tvec = prb.time_grid(tau,x,u);

% Simulate on phyiscal time grid
[~,x2,~] = disc.simulate_dyn(xbar(:,1),{tvec,[u(1:3,:);ones(1,prb.Kfine)]},@(t,x,u) prb.dyn_func(t,x,u),[0,tvec(end)],prb.Kfine,prb.disc);

m = x(1,:);
rI = x(2:4,:);
vI = x(5:7,:);
qBI = x(8:11,:);
omgB = x(12:14,:);

fprintf('\nFinal position error: %.3f\nFinal velocity error: %.3f\n',norm(rI(:,end)-prb.rIK),norm(vI(:,end)-prb.vIK));

save('recent_solution','m','rI','vI','qBI','omgB','tvec','tau','u','x','x2','prb',...
                       'xbar','ubar','tvecbar');

plot_solution;