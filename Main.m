clc;
clear all;
 
global n % number of simulation points
global d % number of decision variables
global d_1 % number of uncertainty variables

% Utility options
plot_kriging = 0
cross_validation_ = 0
%

d=2;  % number of decision variables
d_1=2; % number of noise variables 
n=21; % number of simulation/evaluation points 

global x_mat %simulation points n*d
global w_mean
global e_mat

% standardization of experiment points
% x_mat=(x_mat-repmat(mid_x,n,1))./repmat(rang_x/2,n,1);
% e_mat=(e_mat-repmat(mid_e,n,1))./repmat(rang_e/2,n,1);

L=LHS(-ones(1,d),ones(1,d),-ones(1,d_1),ones(1,d_1),n);
x_mat=L.x;
e_mat=L.e;

% n=size(x_mat,1);  % number of simulation (experiment) points
% d=size(x_mat,2);  % d is the dimension of the decision variables
% d_1=size(e_mat,2); % d_1 is the dimension of the random variables

w_mean=zeros(n,1); % sim outputs mean n*1 (1 output)
%w_var=zeros(n,1); % sim outputs variance n*1 (1 output)


for i=1:size(x_mat,1)
    sim=simulation(x_mat(i,:),e_mat(i,:));
    w_mean(i)=sim;
end

%% variogram
% this code is for n*1
%  sum=zeros(n,1);
%  vario=zeros(n,1);
%  for i=1:n
%      
%      for
%      if (Q_mat(
%      while(Q_mat(j)<=high-i*rang/(n))
%          sum(i)= sum(i)+ (w_mean(j+i)-w_mean(j))^2;
%          j=j+1;
%      end
%      vario(i)=1/(2*j)*sum(i);
%  end
% disp(vario);
% hs=1:1:(n-1);
% plot(hs*rang/(n-1),vario(hs));
%  

%% calling kriging
% example
% x0=[-.4 .3]; %prediction point example 1*d
% e0=[-.2 .3]; %1*d_1
% krig=kriging( (x0-mid_x)./(rang_x/2) , (e0-mid_e)./(rang_e/2))
% sim=simulation(x0,e0)
% [double(krig-sim) double(krig/sim) double(krig) sim]


if plot_kriging==1
	%Kriging %Plot up to d=2
	step=.4; % step of predicting for kriging plot size=1*d Note:steps should creat same number of each decision variable
	step_e=.4;
	x=[];e=[];
	x(:,1)=-.95:step:.95;
	x(:,2)=-.95:step:.95;
	e(:,1)=-.95:step_e:.95;
	e(:,2)=-.95:step_e:.95;
	y=[];
	z=[];
	 for i=1:size(x,1)
		 for j=1:size(e,1)
			krig=kriging( x(i,:), e(j,:));
			y(i,j)=krig;
			sim_p=simulation( x(i,:)*5, e(j,:)*5);
			z(i,j)=sim_p;
		 end
	 end
	surface(1:1:size(x,1),1:1:size(x,1),y);
	hold on
	surface(1:1:size(x,1),1:1:size(x,1),z);
	xlabel('x');
	ylabel('e');
	grid on;
end

if cross_validation_==1
	% Cross Validation

	n=n-1;
	save_x=x_mat;
	save_e=e_mat;
	save_w=w_mean;
	cross=zeros(n-1,2);% [simulation_output kriging_removed_output] oldn-2 removes
	for i=2:n %removing experiment points excluding first and last one (new n is used)
		
		x_mat=save_x; %restoring Q_mat
		e_mat=save_e; %restoring e_mat
		w_mean=save_w; %restoring e_mat
		
		removed_Q=x_mat(i,:);
		removed_e=e_mat(i,:);
		removed_w=w_mean(i,:);
		
		x_mat=x_mat([1:i-1,i+1:end],:);
		e_mat=e_mat([1:i-1,i+1:end],:);
		w_mean=w_mean([1:i-1,i+1:end],:);
		
		krig=kriging( removed_Q ,removed_e ); %fitting kriging on new exp points
		cross(i-1,:)=horzcat(removed_w,krig);
	end

	x_mat=save_x; %restoring Q_mat
	e_mat=save_e; %restoring e_mat
	w_mean=save_w; %restoring e_mat
	n=n+1;
	% Opt 
	non Robust
	x2_start=32000;
	[x,fval]=fmincon( matlabFunction(krig.pred) , x2_start);
	disp([x,fval]);
end


%% PSO
% x0=zeros(1,d); %decision variable vector (d=dimension)
% krig=kriging(x0);
problem.CostFunction = @(x0,e) kriging(x0,e);  % Cost Function
problem.nVar = d;      % Number of Decision Variables
problem.mVar=d_1;      % Number of Noise Variables
problem.VarMin = -ones(1,d);   % Lower Bound of Stadard Decision Variables
problem.VarMax = ones(1,d);   % Upper Bound of Standard Decision Variables
problem.NoiseMax=ones(1,d_1);  % Upper Bound of Stadard Noise Variables
problem.NoiseMin=-ones(1,d_1); % Lower Bound of Stadard Noise Variables
%% Parameters of PSO

% Constriction Coefficients
kappa = 1;
phi1 = 2.05;
phi2 = 2.05;
phi = phi1 + phi2;
chi = 2*kappa/abs(2-phi-sqrt(phi^2-4*phi));

params.PsoMaxIt = 5;       % Maximum Number of Iterations for PSO
params.AlgoMaxIt=45;        % Maximum Number of Iterations for algorithm
params.nPop = 20;           % Population Size (Swarm Size)
params.w = chi;             % Intertia Coefficient
params.wdamp = 1;           % Damping Ratio of Inertia Coefficient
params.c1 = chi*phi1;       % Personal Acceleration Coefficient
params.c2 = chi*phi2;       % Social Acceleration Coefficient
params.ShowIterInfo = true; % Flag for Showing Iteration Informatin

%% Calling PSO

out = PSO(problem, params);

BestSol = out.BestSol;
BestCosts = out.BestCosts;
lhs=out.lhs;

%% Results

figure;
% plot(BestCosts, 'LineWidth', 2);
semilogy(BestCosts, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;


