% Jungho Kim, Junho Song
% Quantile surrogates and sensitivity by adaptive Gaussian process
% for efficient reliability-based design optimization,
% Mechanical systems and signal processing (2021) Vol.161, 107962.
% https://doi.org/10.1016/j.ymssp.2021.107962

clear; close all; tic;

addpath('NIGP/');
addpath('NIGP/util/');
addpath('NIGP/tprod/');

%% Set initial parameters

global nd nc bt stdx Probfun dist_type ub lb model1 model2 history iter...
    comp_id f_opt_sv hyp_optim_o meanfunc covfunc likfunc d_doe Objval

nd = 2;              % number of design parameters
nc = 2;              % number of probabilistic constraints
bt = 2.0;            % target reliability index
dist_type = [1,1];   % distribution type
stdx = [0.1 0.1];    % standard deviation of random variables
lb = [0, 0];         % lower bound
ub = [3.7, 4];       % upper bound
x0 = [1.2, 2.5];     % starting point
rng(11)              % for reproducibility
N_iniDoE = 10;       % number of initial DoE
Opt_size = 2;        % optimization step size
AL_size = 2;         % adaptive learning step size
max_iter = 1e3;      % maximum iteration
eps_f = 1e-3; eps_p = 1e-3;  % tolerance

Costfun = @obj_func;
Probfun = @prob_func;  
BOafun = @Bayopt_func; 
QuanCfun = @prob_func_fmc;

%% Specify the parameters of GP model

meanfunc = [];
covfunc = @covSEard; ell1 = 0.8; ell2 = 0.8; sf = 1; hyp.cov = log([ell1 ell2 sf]);
likfunc = @likGauss; sn = 1e-3; hyp.lik = log(sn);

%% Construct the initial DoE

Coded_value_d = lhsdesign(N_iniDoE,nd);
DoE_bound = [lb+0.5; ub-0.5]';
d_doe = zeros(size(Coded_value_d));
for i = 1:size(Coded_value_d,2)   % Convert coded values to real-world units
    zmax = max(Coded_value_d(:,i));
    zmin = min(Coded_value_d(:,i));
    d_doe(:,i) = interp1([zmin zmax],DoE_bound(i,:),Coded_value_d(:,i));
end
x_doe = d_doe + stdx.*randn(N_iniDoE,nd); % initial DoEs (x)

Objval = Costfun(d_doe);  % objective (cost) function
for i=1:nc  % performance function
    G_val(:,i) = Probfun(x_doe,i) + sn*randn(N_iniDoE,1);
end

%% RBDO by QS2-AGP

iter = 1;
% Set up shared variables with outfun
history.x = [];
history.fval = [];
comp_id_sv=[];
d_doe1 = d_doe; x_doe1 = x_doe;
d_doe2 = d_doe; x_doe2 = x_doe;
G_val1 = G_val(:,1);
G_val2 = G_val(:,2);
fD_stop = 1e4; dD_stop = 1e4;
while (iter <= max_iter)
    %%% adaptive training of quantile surrogates
    hyp_optim_o = minimize(hyp, @gp, -50, @infGaussLik, meanfunc, covfunc, likfunc, d_doe, Objval);
    for k=1:AL_size
        hyp_optim_g1 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, x_doe1, G_val1);
%         hyp_optim_g2 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, x_doe2, G_val2);
%         hyp_optim_g1 = hyp;
        hyp_optim_g2 = hyp;
        model1 = Modelstruct(hyp_optim_g1,x_doe1,G_val1);
        model2 = Modelstruct(hyp_optim_g2,x_doe2,G_val2);
        
        % Find the best design point by quantile surrogates
        options_BO = optimoptions('ga');    % 'Display','iter'
        [d_star, f_star, exitflag_BO, output_BO] = ga(BOafun, nd ,[],[],[],[],lb,ub,[],options_BO);
        BOafun(d_star); % to save comp_id (global variable)
        comp_id_sv = [comp_id_sv; comp_id];
        
        % location of the performance function evaluation
        N_mc = 1e4;
        MCsample_Gp = mvnrnd_rbdo(d_star,N_mc);
        if comp_id==1
            [mu_g_mc,s2_g_mc] = gp(hyp_optim_g1, @infGaussLik, meanfunc, covfunc, likfunc, x_doe1, G_val1, MCsample_Gp);
        else % comp_id==2
            [mu_g_mc,s2_g_mc] = gp(hyp_optim_g2, @infGaussLik, meanfunc, covfunc, likfunc, x_doe2, G_val2, MCsample_Gp);
        end
        U_val_x = gp_U_eval(mu_g_mc, s2_g_mc);
        [U_min_x, i_min_x] = min(U_val_x);
        x_star = MCsample_Gp(i_min_x,:);  % next eval point
        if comp_id==1
            G_star = Probfun(x_star,1) + exp(hyp_optim_g1.lik).*randn(1,1);
        else % comp_id==2
            G_star = Probfun(x_star,2) + exp(hyp_optim_g2.lik).*randn(1,1);
        end

        % Enrich the DoE
        if comp_id==1
            x_doe1 = vertcat(x_doe1,x_star);
            G_val1 = vertcat(G_val1, G_star);
        elseif comp_id==2
            x_doe2 = vertcat(x_doe2,x_star);
            G_val2 = vertcat(G_val2, G_star);
        end
        x_doe = vertcat(x_doe,x_star); % for save
    end
    
    %%% design optimization with parameter sensitivity 
    if iter==1
        d_opt_sv = x0;
        f_opt_sv = Costfun(x0);
        xp = x0;
        [c,~,GC] = QuanCfun(xp);
    else
        xp = d_opt;
    end
    options = optimoptions('fmincon','Display','iter','SpecifyObjectiveGradient',true,...
        'SpecifyConstraintGradient',true,'OutputFcn',@outfun,'MaxIterations',Opt_size);
    [d_opt, f_opt, exitflag, output] = fmincon(Costfun, xp ,[],[],[],[],lb,ub,QuanCfun,options);
    [c,~,GC] = QuanCfun(d_opt);
    
    % current optimum
    d_opt_sv = [d_opt_sv; d_opt]; f_opt_sv = [f_opt_sv; f_opt];
    d_doe = vertcat(d_doe,d_opt); Objval = vertcat(Objval,f_opt);
    fD_stop = abs((f_opt_sv(end-1,:)-f_opt_sv(end,:))./f_opt_sv(end,:));
    dD_stop = abs((d_opt_sv(end-1,:)-d_opt_sv(end,:))./d_opt_sv(end,:));
    if fD_stop <= eps_f && max(dD_stop) <= eps_p
        break
    end
    fprintf('Iteration : %d\n\n',iter); iter = iter+1;
end
iter = iter - 1; toc;

%%% Plot results
grid_interv = 0.05;  % grid interval
[xs1, xs2] = meshgrid(lb(1):grid_interv:ub(1),lb(2):grid_interv:ub(2));
d_plot = [xs1(:) xs2(:)];

Plot_precalc;
Plot_surface;
