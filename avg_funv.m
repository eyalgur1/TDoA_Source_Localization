% Generates and displays a summarizing table that compares the average
% output function value (meaning the value funv(end))


%% Set Parameters
clear
Sigma = logspace(-4, -1, 10)';  % the logspace used in the experiments
sigma_for_table = [-4, -2, -1];  % only the required values of Sigma
N = 15;  % number of sensors in the experiments (4 or 15)
methods = {'TNAM', 'FP', 'SLS', 'SCWLS'};  % if N=15 then SDP cannot be used
s = [2,5];  % the required values of s for T-NAM out of the values used in the experiments
num_nets = 200;  % number of nets used in the experiments
geo = "near";  % array geometry that can be "near", "far" or "circ"

% Values for table (set automatically)
funv_final_table = zeros(length(s) + length(methods) - 1, length(sigma_for_table));  % allocate table
%load(cd+"\output\"+string(N)+geo+"\output_TDOA_"+string(N)+"sen_"+geo+"_TNAM_FP_SOLVIT.mat"); close all  % load the dataset of the iterative methods
TNAM_exists = sum(contains(methods, 'TNAM'));
num_of_methods = TNAM_exists*(length(methods) - 1 + length(s)) + not(TNAM_exists)*length(methods);
names_to_legend = cell(num_of_methods, 1);


%% Create the Table
for j = 1:length(sigma_for_table)
    sigma = find(Sigma == 10^(sigma_for_table(j)));  % loop only over required values of sigma
    current_method = 0;

    for m = 1:length(methods)
        method = methods{m};

        switch method

            case 'TNAM'
                load(cd+"\output\"+string(N)+geo+"\output_TDOA_"+string(N)+"sen_"+geo+"_TNAM_FP_SOLVIT.mat"); close all
                for ss = 1:length(s)
                    current_method = current_method + 1;
                    ssF = s(ss);
                    names_to_legend{current_method} = ['TNAM, s=',num2str(ssF)];
                    for k = 1:num_nets
                        funv_final_table(current_method, j) = funv_final_table(current_method, j) + (1/num_nets)*output.(['net', num2str(k)]).(['sigma', num2str(sigma)]).(method).(['sF',num2str(ssF)]).avg_ML_val(end);
                    end
                end

            case  {'FP', 'SOLVIT'}
                current_method = current_method + 1;
                names_to_legend{current_method} = method;
                for k = 1:num_nets
                    funv_final_table(current_method, j) = funv_final_table(current_method, j) + (1/num_nets)*output.(['net', num2str(k)]).(['sigma', num2str(sigma)]).(method).avg_ML_val(end);
                end

            case  {'SLS', 'SCWLS'}
                load(cd+"\output\"+string(N)+geo+"\output_TDOA_"+string(N)+"sen_"+geo+"_WLS_SCWLS.mat"); close all
                current_method = current_method + 1;
                names_to_legend{current_method} = method;
                for k = 1:num_nets
                    funv_final_table(current_method, j) = funv_final_table(current_method, j) + (1/num_nets)*output.(['net', num2str(k)]).(['sigma', num2str(sigma)]).(method).avg_ML_val(end);
                end
                load(cd+"\output\"+string(N)+geo+"\output_TDOA_"+string(N)+"sen_"+geo+"_TNAM_FP_SOLVIT.mat"); close all

            case 'SDP'
                if N == 4  % SDP exists only for N=4
                    load(cd+"\output\"+string(N)+geo+"\output_TDOA_"+string(N)+"sen_"+geo+"_SDP.mat"); close all  % load the SDP dataset
                    current_method = current_method + 1;
                    names_to_legend{current_method} = method;
                    for k = 1:num_nets
                        funv_final_table(current_method, j) = funv_final_table(current_method, j) + (1/num_nets)*output.(['net', num2str(k)]).(['sigma', num2str(sigma)]).(method).avg_ML_val(end);
                    end
                end
        end
    end
end


%% Generate and Display the Table
avg_out_funv = table(string(names_to_legend), array2table(funv_final_table,'VariableNames',"sigma=10^"+string(sigma_for_table)), 'VariableNames', ["Method", "Avg Output Function Value, N="+string(N)+", "+geo]);
display(avg_out_funv)