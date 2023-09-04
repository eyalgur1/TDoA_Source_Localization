% Generates and displays a summarizing table that compares the average
% obtained RMSE value for each noise levels


%% Set Parameters
clear
Sigma = logspace(-4, -1, 10)';  % the logspace used in the experiments
methods = {'TNAM', 'FP', 'SOLVIT'};  % the correct order is: TNAM -> FP -> SOLVIT -> others (the others can also be without the former three)
s = [1,2,3,4,5];  % the required values of s for T-NAM out of the values used in the experiments
geo = "far";  % array geometry that can be "near", "far" or "circ"
N = 15;

% Values for table (set automatically)
RMSE_table = zeros(length(s) + length(methods), length(Sigma));  % allocate table
TNAM_exists = sum(contains(methods, 'TNAM'));
num_of_methods = TNAM_exists*(length(methods) - 1 + length(s)) + not(TNAM_exists)*length(methods);
names_to_legend = cell(num_of_methods, 1);


%% Create the Table
load(cd+"\output\"+string(N)+geo+"\output_TDOA_"+string(N)+"sen_"+geo+"_TNAM_FP_SOLVIT.mat"); close all  % load the dataset of the iterative methods

for j = 1:length(Sigma)
    sigma = Sigma(j);
    current_method = 0;

    for m = 1:length(methods)
        method = methods{m};

        switch method

            case 'TNAM'
                for ss = 1:length(s)
                    current_method = current_method + 1;
                    ssF = s(ss);
                    names_to_legend{current_method} = ['TNAM, s=',num2str(ssF)];
                    RMSE_table(current_method, j) = output.stat.(['sigma', num2str(j)]).(method).(['sF', num2str(ssF)]).avg_RMSE;
                end
                RMSE_table(end, j) = output.stat.(['sigma', num2str(j)]).avg_CRLB;

            case  {'FP', 'SOLVIT'}
                current_method = current_method + 1;
                names_to_legend{current_method} = method;
                RMSE_table(current_method, j) = output.stat.(['sigma', num2str(j)]).(method).avg_RMSE;

            case {'SLS', 'SCWLS'}
                load(cd+"\output\"+string(N)+geo+"\output_TDOA_"+string(N)+"sen_"+geo+"_WLS_SCWLS.mat"); close all
                current_method = current_method + 1;
                names_to_legend{current_method} = method;
                RMSE_table(current_method, j) = output.stat.(['sigma', num2str(j)]).(method).avg_RMSE;

            case {'SDP'}
                current_method = current_method + 1;
                names_to_legend{current_method} = method;
                if N == 4
                    load(cd+"\output\"+string(N)+geo+"\output_TDOA_"+string(N)+"sen_"+geo+"_SDP.mat"); close all
                    RMSE_table(current_method, j) = output.stat.(['sigma', num2str(j)]).(method).avg_RMSE;
                else
                    RMSE_table(current_method, j) = -1;  % indication that the method did not run for N=15
                end
            load(cd+"\output\"+string(N)+geo+"\output_TDOA_"+string(N)+"sen_"+geo+"_TNAM_FP_SOLVIT.mat"); close all  % load the dataset of the iterative methods
        end
    end
end


%% Generate and Display the Table
names_to_legend{end + 1} = 'CRLB';
T = array2table(RMSE_table,'VariableNames',"sigma=10^"+string(log10(Sigma)));
T_rmse = [table(string(names_to_legend)), T];
disp("Estimated RMSE and CRLB, N="+string(N)+", "+geo)
display(T_rmse)