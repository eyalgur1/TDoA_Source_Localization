% Generates and displays a summarizing table that compares the average 
% obtained function value of the iterative methods at required time stamps


%% Set Parameters
clear
Sigma = logspace(-4, -1, 10)';  % the logspace used in the experiments
sigma_for_table = [-4, -2, -1];  % only the required values of Sigma
N = 4;  % number of sensors in the experiments (4 or 15)
methods = {'TNAM', 'FP', 'SOLVIT'};  % only iterative methods
s = [1,2,3,4,5];  % the required values of s for T-NAM out of the values used in the experiments
num_nets = 200;  % number of nets used in the experiments
geo = "near";  % array geometry that can be "near", "far" or "circ"
intervals = 0:0.00005:0.04;  % the time interval used in the experiments
time_stamp = [0.003, 0.006];  % the required time stamps for funv calculations in a row vector (the stamps must be a member of intervals)

% Values for table (set automatically) 
time_stamp_ind = find(ismember(intervals, time_stamp));
sts = size(time_stamp, 2);  % number of time stamps
funv_time_table = zeros(length(s) + length(methods) - 1, length(sigma_for_table)*sts);  % allocate table
load(cd+"\output\"+string(N)+geo+"\output_TDOA_"+string(N)+"sen_"+geo+"_TNAM_FP_SOLVIT.mat"); close all  % load the dataset of the iterative methods
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
                for ss = 1:length(s)
                    current_method = current_method + 1;
                    ssF = s(ss);
                    names_to_legend{current_method} = ['TNAM, s=',num2str(ssF)];
                    for k = 1:num_nets
                        funv_time_table(current_method, sts*j - sts + 1: sts*j) = funv_time_table(current_method, sts*j - sts + 1: sts*j) + (1/num_nets)*output.(['net', num2str(k)]).(['sigma', num2str(sigma)]).(method).(['sF',num2str(ssF)]).avg_ML_cum(time_stamp_ind);
                    end
                end

            case  {'FP', 'SOLVIT'}
                current_method = current_method + 1;
                names_to_legend{current_method} = method;
                for k = 1:num_nets
                    funv_time_table(current_method, sts*j - sts + 1: sts*j) = funv_time_table(current_method, sts*j - sts + 1: sts*j) + (1/num_nets)*output.(['net', num2str(k)]).(['sigma', num2str(sigma)]).(method).avg_ML_cum(time_stamp_ind);
                end
        end
    end
end


%% Generate and Display the Table
in_table_col = "sigma=10^"+string(kron(sigma_for_table, ones(1, sts)))+", time="+string(repmat(time_stamp, 1, length(sigma_for_table)));
T = table(string(names_to_legend),'VariableNames',"Method");
avg_stamp_funv = [T array2table(funv_time_table,'VariableNames',in_table_col)];
disp("Avg Function Value for Selected Time Stamps, N="+string(N)+", "+geo)
display(avg_stamp_funv)