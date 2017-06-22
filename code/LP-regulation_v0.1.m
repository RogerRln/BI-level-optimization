% The below Optimization Problem is designed with the aim of minimizing the Error between the Predicted
% and Observed Growth Rates while using the less Gene Expression regulation.
%  
%
% OPTIMIZATION
%	Minimize Sum of (growth_predicted_i - growth_observed_i)^2
%	   
%		s.t.  S_i*v_i = 0
%		-(1-indicator)*10^6 <= v_i - v_baseline * GeX <= (1-indicator)*10^6
%
%                i = {1,2,3,..,10}
%		 indicator = {0,1}
%

% LP problem Structure
%  A      LHS matrix
%  b      RHS vector
%  c      Objective coeff vector
%  lb     Lower bound vector
%  ub     Upper bound vector
%  osense Objective sense (-1 max, +1 min)
%  csense Constraint senses, a string containting the constraint sense for
%         each row in A ('E', equality, 'G' greater than, 'L' less than).

% I display the the below example of how to concatenate matrices
% because this is what we do during the first part of the LP formulation
% M = [matrix_1 matrix_2] -> Concatenate Horizontally
% M = [matrix_1; matrix_2] -> Concatenate Vertically

model = readCbModel('gb-2009-10-6-r69-s4', 10000);

%LP.A

% Filling the top left section of the LSH matrix with the 5 Stoichiometric Matrices in the diagonal

zero = zeros(size(model.S)); % Matrix of zeros with dimensions (1383, 1685)

m_cond1 = [model.S zero zero zero zero]; 
m_cond2 = [zero model.S zero zero zero];
m_cond3 = [zero zero model.S zero zero];
m_cond4 = [zero zero zero model.S zero];
m_cond5 = [zero zero zero zero model.S];
%m_cond6 = [zero zero zero zero zero model.S zero zero zero zero];
%m_cond7 = [zero zero zero zero zero zero model.S zero zero zero];
%m_cond8 = [zero zero zero zero zero zero zero model.S zero zero];
%m_cond9 = [zero zero zero zero zero zero zero zero model.S zero];
%m_cond10 = [zero zero zero zero zero zero zero zero zero model.S];

matrix_conditions = [m_cond1; m_cond2; m_cond3; m_cond4; m_cond5];

% Filling the top right section of the LSH matrix zeros

matrix_conditions = [matrix_conditions zeros(size(matrix_conditions))];


% filling the middle left section of the LSH matrix (Identity matrices)

I = eye(size(model.S,2)); % Identity Matrix with dimensions (1685, 1685)

zero = zeros(size(I)); % Matrix full of zeros with dimensions (1685, 1685)

m_identity1 = [I zero zero zero zero]; 
m_identity2 = [zero I zero zero zero];
m_identity3 = [zero zero I zero zero];
m_identity4 = [zero zero zero I zero];
m_identity5 = [zero zero zero zero I];
%identity6 = [zero zero zero zero zero I zero zero zero zero];
%identity7 = [zero zero zero zero zero zero I zero zero zero];
%identity8 = [zero zero zero zero zero zero zero I zero zero];
%identity9 = [zero zero zero zero zero zero zero zero I zero];
%identity10 = [zero zero zero zero zero zero zero zero zero I];

identity_matrix = [m_identity1; m_identity2; m_identity3; m_identity4; m_identity5];

% filling the middle right section of the LSH matrix (GeX ratio matrices) 

%Reading the file with GeX ratios for 5 conditions
gex_ratios = 'gex_ratios.txt';

data = load(gex_ratios, '-ascii');
gex_ratio1 = data(:,1); % column 1 with the Gex ratios of condition 1
gex_ratio2 = data(:,2); 
gex_ratio3 = data(:,3);
gex_ratio4 = data(:,4);
gex_ratio5 = data(:,5);

% Reading the file with the Enzymatic rxns identifiers
enz_rxns = textread('enzymatic_rxns.txt', '%s', 'delimiter', '\n');

Identity_1 = I;
Identity_2 = I;
Identity_3 = I;
Identity_4 = I;
Identity_5 = I;
% Changing the One's of the Identity Matrices to GeX ratios of the corresponding conditions
for i = 1:length(enz_rxns)
	j= find(strcmp(model.rxns, enz_rxns(i)) == 1);
	Identity_1(j,j) = -gex_ratio1(i);
	Identity_2(j,j) = -gex_ratio2(i);
	Identity_3(j,j) = -gex_ratio3(i);
	Identity_4(j,j) = -gex_ratio4(i);
	Identity_5(j,j) = -gex_ratio5(i);
end

m_GeX1 = [Identity_1 zero zero zero zero];
m_GeX2 = [zero Identity_2 zero zero zero];
m_GeX3 = [zero zero Identity_3 zero zero];
m_GeX4 = [zero zero zero Identity_4 zero];
m_GeX5 = [zero zero zero zero Identity_5];

GeX_matrix = [m_GeX1; m_GeX2; m_GeX3; m_GeX4; m_GeX5];

% When I was trying to concatenate the below two matrices (with the 10 conditions) I was getting an OUT OF MEMORY Error, that is why I reduced the problem to five conditions.

identity_GeX = [identity_matrix GeX_matrix];

% filling the low section of the LSH matrix (identity_GeX)
% Getting OUT OF MEMORY ERROR
identity_GeX = [identity_GeX;identity_GeX];


%LP.c
% First, we add an artificial metabolite that consists of putting stoichiometric coefficients of 1 in all the alpha columns.
new_matrix(end+1,:) = zeros(size(new_matrix(1,:)));
alfas = 1686:3370;
new_matrix(end,alfas) = +1;

% Second, we add an artificial reaction that "consumes" the above artificial metabolite and mimics the sum of alphas.
new_matrix(:,end + 1) = zeros(size(new_matrix(:,1)));
new_matrix(end, end) = -1;

%Third, as we want to minimize the sum of alphas, we select the above reaction as our new objective function (we later minimize this objective function during the optimization).
c = [model.c;model.c;1];
c(427) = 0;
c(2112) = 0;


%LP.lb
alfas_zeros = zeros(size(model.lb));
LB = [model.lb; alfas_zeros];
LB(end+1) = 0;


%LP.ub
alfas_ones= ones(size(model.lb));
if optimization ~= 1
	alfas_ones(1:1685,1) = 1000;
end
UP = [model.ub; alfas_ones];
UP(end+1) = 10000;


%LP.csense
[nMets,nRxns] = size(model.S);
csense_equality(1:nMets,1) = 'E';
csense_greater(1:nRxns,1) = 'G';
csense_less(1:nRxns,1) = 'L';
csense =[csense_equality; csense_greater; csense_less];
csense(end+1) = 'E';

%Lp.b
b_stoimatrix(1:nMets,1) = 0;
b_upper = modelMedium.ub;
b_lower = modelMedium.lb;
B =[b_stoimatrix; b_lower; b_upper];
B(end+1) = 0;

%%%%%%%%%%%%%%%%%%%%%%%% Optimization 1 %%%%%%%%%%%%%%%%%%%%%%

if optimization == 1
	%LP.vtype   % Reactions are treated as continuous Variables and Alphas as Binary Variables
	[nMets,nRxns] = size(model.S);
	vtype_rxns = repmat('C', nRxns, 1);
	vtype_alphas = repmat('B', nRxns, 1);
	vtype =[vtype_rxns; vtype_alphas];
	vtype(end+1) = 'C';

	number_rxns=[];

	for i = 1:length(gr)
    		% We fix the growth rate
    		LB(427) = gr(i);
    		UP(427) = gr(i);
    		% Minimize the objective function (Sum of Alphas)
    		LP.osense = +1;
    		LP.A = new_matrix;
    		LP.c = c;
    		LP.lb = LB;
    		LP.ub = UP;
    		LP.csense = csense;
    		LP.b = B;
    		LP.vtype = vtype;
    		solution = solveCobraLP(LP);
    		n_rxns = solution.obj;
    		number_rxns = [number_rxns n_rxns];
	end

	rxns = [];
	percentage_change = [];
	vector_distribution_pFBA = [];
	vector_distribution_FBA = [];

	LP =[];


    	s = 'number_rxns_opt1_SMMcutoff.txt';
    	fileID = fopen(s,'w');
    	for row = 1:length(number_rxns)
        	fprintf(fileID,  '%.2f\n', number_rxns(row));
    	end
    	fclose(fileID);

end
