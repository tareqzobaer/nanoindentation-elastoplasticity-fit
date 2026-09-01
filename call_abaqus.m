test_file_names = {'sp_ind_sample1_1.xlsx','sp_ind_sample1_2.xlsx','sp_ind_sample2.xlsx'...
    'sp_ind_sample3.xlsx','sp_ind_sample5_1.xlsx','sp_ind_sample5_2.xlsx','sp_ind_sample6_1.xlsx','sp_ind_sample7.xlsx'};

test_files = 1:length(test_file_names);

for j=test_files
    tag = test_file_names{j}(8:end-4);
    test_file_name = strcat('ind_se_data\',test_file_names{j});
    ind_file_name = strcat('ind_se_data\',test_file_names{j}(1:7),'elastic_',tag,'csv');

    indst_data = readmatrix(ind_file_name);

    datafile = strcat('ind_se_data\p_',tag,'mat');
    if ~isfile(datafile)
        p_fit(tag);
    end
    load(datafile,'p','n','umax', 'sheet_names');
    test_num = length(sheet_names);


    tests = 1:test_num;

    for i=tests
        sheet = sheet_names{i};
        sheet_num = str2double(sheet(end-2:end));
        nu = 0.3;
        additional_args = {nu, -umax(i), p{i},n};

        ind_data = indst_data(sheet_num,2:end); % [E, y0, e0]

        ind_data(:,1) = [];
        data_variation = [0 1 2.5 5 10];
        n_ind = 2;

        % GA parameters
        nvars = 6; % [y0,e0,y1,e1,y2,e2]
        total_popsize = 50;
        MaxGenerations = 15;
        MaxStallGenerations = 5;

        lb = [ 0.9*ind_data(1), 0.9*ind_data(2), ind_data(1), ind_data(2), ind_data(1), ind_data(2)];   % Lower bounds
        ub = [  1.1*ind_data(1), 1.1*ind_data(2), 300, 3.0/100, 300, 3.0/100];   % Upper bounds

        A = [1 0 -1  0  0  0;  % y0 - y1 < 0
            0 0  1  0 -1  0;  % y1 - y2 < 0
            0 1  0 -1  0  0;  % e1 - e2 < 0
            0 0  0  1  0 -1];  % e0 - e1 < 0
        b= zeros(4, 1);

        [biased_pop] = create_biased_init_pop(ind_data, data_variation, n_ind, lb(3:6), ub(3:6), nvars, A, b);


        popsize = total_popsize - n_ind*length(data_variation);
        population_ul = zeros(popsize, nvars);
        invalid_data = true(popsize, 1);
        length_ivd = nnz(invalid_data);
        counter = 0;
        while length_ivd~=0

            population_ul(invalid_data,:) = (repmat(ub, length_ivd, 1)...
                - repmat(lb, length_ivd, 1)).*lhsdesign(length_ivd,nvars)...
                + repmat(lb, length_ivd, 1);

            for k=find(invalid_data)
                adata = all(  A*population_ul(k,:)'  <b);
                invalid_data(k) = ~(adata); 
            end

            length_ivd = nnz(invalid_data);
            counter = counter+1;
            if counter>10000
                break
            end
        end

        % Create the initial population matrix
        population = [population_ul; biased_pop];
        population = population(randperm(total_popsize),:);
        ranges = ub - lb;

        % Normalize the initial population
        population_normalized = bsxfun(@minus, population, lb) ./ ranges;
        additional_args_nonlincon = {ranges,lb};

        % Run the genetic algorithm optimization
        objective_function = @(x) calculate_mse(x.*ranges + lb, additional_args);
        nonlincon = @(x) trilinear_matmodel(x.*ranges + lb);

        options = optimoptions('ga', 'InitialPopulation', population_normalized, ...
            'Populationsize',total_popsize,'MaxGenerations',MaxGenerations,'MaxStallGenerations',MaxStallGenerations,...
            'UseParallel', false,'Display', 'iter');

        fprintf('%s, Test: %s \n\n', test_file_names{j}(8:end-5),sheet(end-2:end))
        [x_optimal_normalized, fval_optimal] = ga(objective_function, nvars, A, b, [], [], zeros(1, nvars), ones(1, nvars), [], options);
        
        % Denormalize the optimal solution
        x_optimal = x_optimal_normalized .* ranges + lb;

        y0_opt = x_optimal(1);
        e0_opt = x_optimal(2);
        y1_opt = x_optimal(3);
        e1_opt = x_optimal(4);
        y2_opt = x_optimal(5);
        e2_opt = x_optimal(6);


        E0 =  (y0_opt-0)/(e0_opt-0)/1e3;

        % Save results
        save_results = strcat('ga_results\',test_file_names{j}(8:end-5),'test',sheet(end-2:end),'.txt');
        fileID = fopen(save_results, 'w');
        fprintf(fileID,'%s',sprintf('File: %s, Test: %s \n', test_file_names{j}(8:end-5),sheet(end-2:end)));
        fprintf(fileID,'%s',newline);
        fprintf(fileID,'%s',sprintf('------------------------------------------------\n'));
        fprintf(fileID,'%s',sprintf('           Linear elasticity data           \n'));
        fprintf(fileID,'%s',sprintf('------------------------------------------------\n'));
        fprintf(fileID,'%s',sprintf('E0 = %0.5f GPa \n',E0));
        fprintf(fileID,'%s',sprintf('y_0 = %0.5f MPa \n',y0_opt));
        fprintf(fileID,'%s',sprintf('e_0 = %0.5f \n',e0_opt));
        fprintf(fileID,'%s',newline);
        fprintf(fileID,'%s',sprintf('nu = %0.5f \n',nu));
        fprintf(fileID,'%s',newline);
        fprintf(fileID,'%s',sprintf('------------------------------------------------\n'));
        fprintf(fileID,'%s',sprintf('           Isotropic plasticity data           \n'));
        fprintf(fileID,'%s',sprintf('------------------------------------------------\n'));
        fprintf(fileID,'%s',sprintf('y_1 = %0.5f MPa \n',y1_opt));
        fprintf(fileID,'%s',sprintf('e_1 = %0.5f \n',e1_opt));
        fprintf(fileID,'%s',sprintf('y_2 = %0.5f MPa \n',y2_opt));
        fprintf(fileID,'%s',sprintf('e_2 = %0.5f \n',e2_opt));
        fprintf(fileID,'%s',newline);
        fprintf(fileID,'%s',sprintf('------------------------------------------------\n'));
        fprintf(fileID,'%s',sprintf('           Mean of Squared Errors          \n'));
        fprintf(fileID,'%s',sprintf('------------------------------------------------\n'));
        fprintf(fileID,'%s',sprintf('MSE = %0.5f \n',fval_optimal));
        fclose(fileID);

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mse] = calculate_mse(x, additional_args)

% {nu, umax, p{i},n};
nu = additional_args{1};
umax = additional_args{2};
p = additional_args{3};
n = additional_args{4};

% Define the paths
baseFile = 'base\abaqus_test9.inp';

jobName = strcat('job_gen_cand_',num2str(randi(10000)+randi(10000)));
pythonScript_base = 'base\get_data_at_RP.py';

newjobinp = edit_inp(baseFile,jobName,x,nu,umax);
pythonScript = edit_py(pythonScript_base,jobName);


% Run Abaqus job
% system(['abaqus job=', jobName, ' input=', newjobinp]);
system(['abaqus job=', jobName, ' input=', newjobinp, ' cpus=6']);

% check log file to see if simulation completed
log_file = strcat(jobName,'.log');
analysis_status = 0;

while ~analysis_status
    [analysis_status]=analysis_state(log_file);
    pause(3)
end
if analysis_status==2 % error encountered
    mse = 1e10;
    delete_files(jobName, analysis_status)
    return
end

% Extract data from odb to csv
[status, cmdout] = system(['abaqus python ', pythonScript]); %#ok<ASGLU>
%
% if status~=0
%     mse = 1e10;
%     return
% end
ph_data = readmatrix(strcat(jobName,'.csv'));
simulation_h = ph_data(:,3).*-1; % um
simulation_P = ph_data(:,4).*-1; % mN

q = polyfit(simulation_h,simulation_P,n);
test_points = 1000;
% figure(1)
% plot(simulation_h,simulation_P,'LineWidth',2)
% hold on
% scatter(simulation_h,simulation_P)
% hold off
% test_h = 0:0.1:simulation_h(end);
test_h = linspace(0,simulation_h(end),test_points);

experimental_P = polyval(p,test_h);
fea_P = polyval(q,test_h);
sse = sum((experimental_P - fea_P).^2);
mse = sse/test_points;
delete_files(jobName, analysis_status)
end

function [c,ceq] = trilinear_matmodel(x)

y0 = x(1);
e0 = x(2);
y1 = x(3);
e1 = x(4);
y2 = x(5);
e2 = x(6);

E0 =  (y0-0)/(e0-0);
Ep1 = (y1-y0)/(e1-e0);
Ep2 = (y2-y1)/(e2-e1);

c(1) = Ep1 - E0;
c(2) = Ep2 - E0;
ceq=[];

% Check for NaN or Inf
if any(isnan(c) | isinf(c))
    % warning('Constraint function produced NaN or Inf');
    % disp(x)
    c(isnan(c) | isinf(c)) = 1e10;  % Replace with a large number, or handle as appropriate
end

end


function [out_pop] = create_biased_init_pop(ind_data, data_variation, n_ind, lb_t12, ub_t12, nvars, A, b)

data_var_l = 1 - data_variation./100;
data_var_u = 1 + data_variation./100;

ind_pop = zeros(n_ind*length(data_variation),length(ind_data));
pop_idx = reshape(1:n_ind*length(data_variation)',n_ind,[])';

for i=1:length(data_variation)
    ind_pop(pop_idx(i,:),:) = repmat(data_var_l(i).*ind_data, n_ind, 1) + ...
        repmat(data_var_u(i).*ind_data - data_var_l(i).*ind_data, n_ind, 1) .* lhsdesign(n_ind, length(ind_data));
end

invalid_data = true(n_ind*length(data_variation), 1);
y_pop = zeros(n_ind*length(data_variation), nvars-length(ind_data));


out_pop = zeros(n_ind*length(data_variation), nvars);
length_ivd = nnz(invalid_data);
counter = 0;
while length_ivd~=0

    y_pop(invalid_data,:) = (repmat(ub_t12, length_ivd, 1)...
        - repmat(lb_t12, length_ivd, 1)).*rand(length_ivd,nvars-length(ind_data))...
        + repmat(lb_t12, length_ivd, 1);

    out_pop(invalid_data,:)  = [ind_pop(invalid_data,:) y_pop(invalid_data, :)];
    for k=find(invalid_data)
        adata = all(  A*out_pop(k,:)'  <b);
        invalid_data(k) = ~(adata);
    end

    length_ivd = nnz(invalid_data);

    if counter>50000
        break
    end
    counter = counter+1;
end
end




function [destinationFile]=edit_inp(inpFile,jobName,x,nu,umax)

% unpack variables
y0 = x(1)/1e3; % MPa -->GPa
yp1 = x(3)/1e3; % MPa -->GPa
yp2 = x(5)/1e3; % MPa -->GPa


e0 = x(2);
e1 = x(4);
e2 = x(6);


E0 =  (y0-0)/(e0-0);

% ep0 = e0-e0;
ep1 = e1-e0;
ep2 = e2-e0;

sourceFile = inpFile;
destinationFile = strcat(jobName,'.inp');
[~] = copyfile(sourceFile, destinationFile);

% Specify the line number you want to modify
linesToModify = [5366 5368 5369 5370 5404]; % test9
modified_prop = {sprintf('%s%s%s%s%s%s',num2str(E0),', ', num2str(nu)),...
    sprintf('%s%s%s%s',num2str(y0),', ', num2str(0),'.'),...
    sprintf('%s%s%s',num2str(yp1),', ', num2str(ep1)),...
    sprintf('%s%s%s',num2str(yp2),', ', num2str(ep2)),...
    sprintf('%s%s','RP_set, 2, 2, ', num2str(umax))};

% Open the file for reading
fid = fopen(destinationFile, 'r');
if fid == -1
    error('Base input file could not be opened for reading');
end

% Read all lines into a cell array
lines = {};
lineIndex = 1;
while ~feof(fid)
    lines{lineIndex} = fgetl(fid); %#ok<AGROW>
    lineIndex = lineIndex + 1;
end

% Close the file after reading
fclose(fid);

for i=1:length(linesToModify)
    lines{linesToModify(i)} = modified_prop{i}; %#ok<AGROW>
end

% Open the file for writing (overwrite mode)
fid = fopen(destinationFile, 'w');
if fid == -1
    error('Job inputfile could not be opened for writing');
end

% Write all lines back to the file
for i = 1:length(lines)
    fprintf(fid, '%s\n', lines{i});
end

% Close the file after writing
fclose(fid);


end

function [analysis_status]=analysis_state(log_file)
fid = fopen(log_file, 'r');
if fid == -1
    % error('Log file could not be opened for reading');
    analysis_status = 2;
    return
end

lines = {};
lineIndex = 1;
while ~feof(fid)
    lines{lineIndex} = fgetl(fid); %#ok<AGROW>
    lineIndex = lineIndex + 1;
end
fclose(fid);

indices = strfind(lines{end}, 'COMPLETED');

if ~isempty(indices)
    % disp('Analysis completed')
    analysis_status = 1;
else
    findice = strfind(lines{end}, 'error');
    if isempty(findice)
        % disp('Analysis running')
        analysis_status = 0;
    else
        % disp('Analysis resulted in error')
        analysis_status = 2;
    end
end

end

function newpythonScript= edit_py(pythonScript,jobName)
sourceFile = pythonScript;
newpythonScript = strcat(jobName,'.py');
[~] = copyfile(sourceFile, newpythonScript);

fid = fopen(newpythonScript, 'r');
if fid == -1
    error('Python file could not be opened for reading');
end
lines = {};
lineIndex = 1;
while ~feof(fid)
    lines{lineIndex} = fgetl(fid); %#ok<AGROW>
    lineIndex = lineIndex + 1;
end
fclose(fid);
lines{6} =  sprintf('%s%s%s%s','odb_path = ','''',jobName,'.odb''');
lines{7} =  sprintf('%s%s%s%s','csv_file = ','''',jobName,'.csv''');

fid = fopen(newpythonScript, 'w');
if fid == -1
    error('Job inputfile could not be opened for writing');
end

% Write all lines back to the file
for i = 1:length(lines)
    fprintf(fid, '%s\n', lines{i});
end

% Close the file after writing
fclose(fid);
end



function delete_files(jobName, analysis_status)
if analysis_status==2
    fileTypes = {'.msg', '.prt', '.sta', '.csv'};
else
    fileTypes = {'.msg', '.com','.dat','.inp','.log','.odb','.prt','.sta','.csv','.py','.sim','.for'};
end
for i = 1:length(fileTypes)
    filename = strcat(jobName, fileTypes{i});
    if isfile(filename)
        delete(filename);
    end
end

end

function p_fit(tag)

test_file_name = strcat('data\sp_',tag,'xlsx');

[~, sheet_names] = xlsfinfo(test_file_name);
results_data = readmatrix(test_file_name, 'Sheet', sheet_names{1}, 'UseExcel', 0);

sheet_names = sheet_names(4:length(sheet_names)-1);
test_num = length(sheet_names);

tests = 1:test_num;
results_data = results_data(tests+1,:);

p = cell(test_num,1);
umax = zeros(test_num,1);
unload_E = zeros(test_num,1);
for i=tests
    sheet = sheet_names{i};
    sheet_num = str2double(sheet(end-2:end));
    unload_E(i) = results_data(sheet_num,2);

    loading_data = readtable(test_file_name, 'Sheet', sheet, 'UseExcel', 0);
    loading_data = loading_data(:,1);
    load_smt = ~cellfun(@isempty,loading_data.Segment);

    segments = find(load_smt);

    load_start= segments(1);
    load_end = segments(2)-1;

    columns = [2 3];

    data = readmatrix(test_file_name, 'Sheet', sheet, 'UseExcel', 0);

    data = data(load_start:load_end,columns);

    nan_data = any(isnan(data),2);
    data(nan_data,:) = [];

    % h_t: total penetration depth (directly measured)
    % P: load (directly measured)


    h_t = data(:,1).*1e-3;
    umax(i) = max(h_t);
    P = data(:,2);
    
    n=12;
    
    [p{i},s] = polyfit(h_t,P,n);
end

save(strcat('ind_se_data\p_',tag,'mat'),'p','n','umax', 'sheet_names')
end

