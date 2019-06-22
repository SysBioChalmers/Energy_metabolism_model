function sol = solveMOMENT(model,objective,osenseStr,prot_cost_info,tot_prot_weight)

[nMets,~] = size(model.S);

% Determine objective
model.c(strcmp(model.rxns,objective)) = 1;

% Construct LP
A = (prot_cost_info)';
b = tot_prot_weight*1000;

Aeq = model.S;
beq = zeros(nMets,1);

lb = model.lb;
ub = model.ub;

%linprog always runs minimization
if osenseStr == 'max'
    f = -model.c;
elseif osenseStr == 'min'
    f = model.c;
end

[x,~,exitflag,~] = linprog(f,A,b,Aeq,beq,lb,ub);

sol = struct();
sol.fluxes = x;
sol.exitflag = exitflag;

if exitflag == 1
	sol.protUsage = (prot_cost_info)' * x / 1000;
end

    