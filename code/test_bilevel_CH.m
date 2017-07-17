% Bi-level Bacillus subt. Problem
% The below Bi-level Problem is designed with the aim of minimizing the difference between the Predicted
% and Observed Growth Rates (of Bacillus subtilis) while Maximizing the biomass production under 10 different conditions.
% We made use of binary variables in order to identify the Gene Regulation that helps to minimize the above difference. 
% 
%  
%
% OPTIMIZATION
%	Minimize ∑(growth_predicted_i - growth_observed_i)^2
%
%		Maximize (growth_predicted_i)
%			s.t.  S_i*v_i = 0
%			-(1-indicator)*10^6 <= v_i - v_baseline * GeX <= (1-indicator)*10^6
%
%                	i = {1,2,3,..,10}
%		 	indicator = {0,1}
%
%	∑(indicators) = fixed number (?)


% Testing the problem in one condition: Growth on CH using SMM as Baseline, we try to identify the set of genes that regulates metabolism
% such that the biomass flux value is close to the Observed growth rate of Bacillus subtilis on CH medium.


% Load Gene Expression Ratio Data set for the CH condition, the Gex Ratios will constrain the bounds of the enzymatically-catalyzed reactions.

gex_ratios = 'gex_ratio_CH.txt';
data = load(gex_ratios, '-ascii');
gex_ratio1 = data(:,1); % column with the Gex ratios


% Reading the file with the enzymatically-catalyzed reactions identifiers
enz_rxns = textread('enzymatic_rxns.txt', '%s', 'delimiter', '\n');


% Blocked transport reactions (amino acids, etc), we need to re-open these reactions.
make_reversible = {'rxn05466', 'rxn05217', 'rxn05305', 'rxn05215', 'rxn00910', 'rxn05244', 'rxn05299', 'rxn00693', 'rxn01816', 'rxn05663', 'rxn05221', 'rxn05669', 'rxn05243'};
make_forward = {'rxn05219', 'rxn05303'};

% Vmax and Vmin results from the FVA of the Baseline condition, These initial lower and upper bounds are multiplied by the Gex Ratios.
vmax = load('Vmax_90_SMM.txt', '-ascii');
vmin = load('Vmin_90_SMM.txt','-ascii');

for i = 1:length(model.rxns)
	if vmin(i) > vmax(i)
		new_vmax = vmin(i);
		vmin(i) = vmax(i);
		vmax(i) = new_vmax;
	end

	upper = max(0, vmax(i));
	lower = min(0, vmin(i));
	vmin(i) = lower;
	vmax(i) = upper;

	if find(strcmp(make_reversible, model.rxns(i)) == 1)
		vmax(i) = 1000;
		vmin(i) = -1000;
	end
	
	if find(strcmp(make_forward, model.rxns(i)) == 1)
		vmax(i) = 1000;
		vmin(i) = 0;
	end
end


% Identify the index of the enzymatically-catalyzed reactions

enzymatic = zeros(length(enz_rxns)-15,1);
change_direction = [make_reversible make_forward];
j=0;
for i = 1:length(enz_rxns)

	ind = find(strcmp(model.rxns, enz_rxns(i)) == 1);

	if sum(strcmp(change_direction, enz_rxns(i))) == 0
		j = j + 1;
		enzymatic(j) = ind;
	end
end


gex_index = zeros(length(enzymatic),1);
for i = 1:length(enzymatic)
	ind = find(strcmp(enz_rxns, model.rxns(enzymatic(i))) ==1);
	gex_index(i) = ind;
end


% Index for reversible reactions
rev = find(model.rev == 1);

% Index for irreversible reactions
irr = find(model.rev == 0);


%%%%%%%%%%%%%%%%%%% Set up and solve the Bi-level optimization problem  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define decision (fluxes) and binary (indicators) variables

v = sdpvar(length(model.rxns),1); % 1685 fluxes
indicators = binvar(length(enzymatic), 1); % 1243 indicators


% Outer objective
OO = (v(427) - 1.2);

% Outer constraint
CO = [sum(indicators) >= 1220];


% Inner objective
OI = -v(427);

% Inner constraint
CI = [model.S*v == 0,
     -10000 <= v(rev) <= 10000
          0 <= v(irr) <= 10000,
     -1e6*indicators  <= v(enzymatic) - vmax(enzymatic).*gex_ratio1(gex_index) <= 1e6*indicators];


% Solve the problem
solvebilevel(CO,OO^2,CI,OI,v)

value(v(427))
