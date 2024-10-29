% % pyOut: This is the similarity matrix that contains the similarity scores between candidate molecules, which are used in optimization.
% A2: This matrix contains collected data with the following columns:
% Column 1-9: Empty, as sharing the experimental data is restricted.
% Column 1-7: Area [%] of individual peaks, representing the contribution of each peak to the total area under the curve.
% Column 8: Carbon number, indicating the number of carbon atoms in the molecular structure.
% Column 9: Retention time, the time taken for a compound to pass through the chromatography column.
% Column 10: Kovats index, a measure for identifying compounds based on their retention time.
% Column 11: Relative retention time between two paraffins, providing a comparative metric for analyzing retention times.
% Column 12: Corresponding SMILES (Simplified Molecular Input Line Entry System) strings for molecular candidates, representing their chemical structure.
% Column 13: Number of candidates for each molecular structure, indicating how many distinct structures correspond to each entry.
% Column 14: Unique identifier for candidate molecules (CID index - PubChem specific), a specific identifier used for tracking molecular candidates.
% Column 15: Deviation of candidates in retention index from the peaks' retention index, showing how much the candidate's retention index varies from expected values.
% Column 16: Number of the peak corresponding to each candidate, linking candidates back to specific peaks in the chromatographic analysis.


clear all
clc
close all

% Global variables for optimization constraints and outputs
global pyOut f2 A Aeq b beq

% Load collected data and similarity matrix
load('CollectedData.mat') % Load A2 containing candidate data
load('similarityMatrix.mat') % Load the similarity matrix into pyOut

a = A2; % Rename the variable for clarity, keeping the original data intact
P =[]; % Initialize the output matrix
id2 = 0; % Index for column assignments in P
there_is_identical_lines = []; % Store lines that have identical identifiers

% Loop through each candidate to create a list of identifiers
for kidl = 1:size(a,1)
    if A2{kidl,13} ~= 1 % If the number of candidates is not 1
        there_is_identical_lines = [there_is_identical_lines; A2{kidl,14}]; % Add CID index
    else
        there_is_identical_lines = [there_is_identical_lines; repmat(A2{kidl,14}, 1, 5)]; % Replicate CID index for candidates
    end
end

keyboard % Pause for debugging

% Identify unique identifiers and duplicates
[u,I,J] = unique(there_is_identical_lines, 'rows', 'first');
hasDuplicates = size(u,1) < size(there_is_identical_lines,1); % Check for duplicates
ixDupRows = setdiff(1:size(there_is_identical_lines,1), I); % Find indices of duplicate rows

% Initialize plus_line for future matrix manipulations
plus_line = zeros(size(pyOut,1), 1); 
for ka = 1: find([A2{:,11}] == 0, 1, 'last') % Loop until the first instance of column 11 being zero
    if any(ixDupRows == ka) % Check if current index is a duplicate
        a{ka, 12}{end + 1} = ['dummy' num2str(ka)]; % Add a dummy entry in the smiles column
        a{ka, 14}(end + 1) = mean(a{ka, 14}) + rand(1); % Add random mean value to CID index
        a{ka, 15}(end + 1) = max([a{:, 15}]); % Add max deviation to the deviation column
        a{ka, 16}(end + 1) = a{ka, 16}(end); % Replicate peak number
        
        id2 = id2 + 5; % Increment id2 for the next column assignment
        P = [P pyOut(:, id2 - 4:id2) plus_line]; % Append similarity matrix columns to P
    else
        if a{ka, 13} ~= 1 % If the number of candidates is not 1
            id2 = id2 + 5; 
            P = [P pyOut(:, id2 - 4:id2)]; % Append similarity columns for multiple candidates
        else
            id2 = id2 + 1; 
            P = [P pyOut(:, id2)]; % Append single similarity column for single candidate
        end
    end
end

% Initialize P2 to organize output results
P2 = [];
id2 = 0; % Reset index for P2
plus_line = zeros(1, size(P,2)); % Initialize a row for padding

for ka = 1: find([A2{:,11}] == 0, 1, 'last') % Loop until the first instance of column 11 being zero
    if any(ixDupRows == ka) % Check if current index is a duplicate
        id2 = id2 + 5; % Increment for duplicates
        P2 = [P2; P(id2 - 4:id2, :); plus_line]; % Append corresponding rows to P2 and add padding
    else
        if a{ka, 13} ~= 1 % If the number of candidates is not 1
            id2 = id2 + 5; 
            P2 = [P2; P(id2 - 4:id2, :)]; % Append multiple candidate rows
        else
            id2 = id2 + 1; 
            P2 = [P2; P(id2, :)]; % Append single candidate row
        end
    end
end

pyOut = P2; % Update the global similarity matrix with the processed values
clearvars P P2 % Clear intermediate variables to save memory

% Extract relevant indices for further processing
num = find([a{:,8}] < 32 & [a{:,8}] > 25, 1, 'last'); % Find paraffin indices
SMLS = [a{1:max(num),12}]; % Extract smiles from the specified range
CIDSS = [a{1:max(num),14}]; % Extract CID indices

% Adjust the similarity matrix size to match CIDSS
pyOut = pyOut(1:length(CIDSS), 1:length(CIDSS));

% Create peak identifiers based on unique values
u = [a{1:max(num),16}]; % Peak identifiers
u_unique = unique(u); % Get unique peak identifiers
for u_k = 1:numel(u_unique) % Map original identifiers to unique values
    u(u == u_unique(u_k)) = u_k; % Replace original peak identifiers
end

keyboard % Pause for debugging

% Create constraints for optimization based on CIDSS and unique peaks
[A,b,Aeq,beq] = create_contraints(CIDSS,u);

% Initialize the objective function
f2 = zeros(1, length([a{1:max(num),15}])); 
idzt = 1; 
for kzt = 1:max(num) % Loop through each candidate
    zt = -[a{kzt,15}]; % Negative deviation values
    min_t = min(zt,[],2); % Minimum deviation
    max_t = max(zt,[],2); % Maximum deviation
    f2(idzt:idzt + length(zt) - 1) = 0.5 * (zt - min_t) ./ (max_t - min_t + eps) + 0.5; % Normalize
    idzt = idzt + length(zt); % Increment index for f2
end
f2 = repmat(f2, 1, size(pyOut,1)); % Replicate f2 to match the size of pyOut

% Define optimization constraints
lb = zeros(1, size(f2,2)); % Lower bounds
ub = ones(1, size(f2,2)); % Upper bounds
intcon = 1:size(f2,2); % Integer constraints
options = optimoptions('intlinprog', 'Display', 'off'); % Optimization options
f = (reshape(-pyOut', [], 1) - f2')'; % Objective function for optimization

clearvars -except actual_name A2 actual_files list kSJ1 kSJ2 Files Files_PTBA names f A b ub lb options intcon Aeq beq CIDSS FN SMLS % Keep necessary variables for further use

% Solve the integer linear programming problem
[x, fval, exitflag, output] = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub, [], options);

clearvars A b Aeq beq lb ub f % Clear variables used for the optimization

% Process the output to identify selected candidates
temp_vector = 1:length(x); % Create a vector for indices
temp_vector(x == 0) = 0; % Set selected indices to zero
temp_matrix = reshape(temp_vector, sqrt(length(x)), []); % Reshape into a square matrix
ID_maindiagonal = find(eye(sqrt(length(x)))); % Find the main diagonal
IDS = temp_matrix(ID_maindiagonal); % Extract diagonal indices
IDS(IDS == 0) = []; % Remove zero entries
IDS_rem = (rem(IDS, sqrt(length(x)))); % Remainder for reshaping
IDS_int = floor(IDS ./ sqrt(length(x))); % Integer division for reshaping
IDS_final = IDS_rem; % Initialize final indices
IDS_final(IDS_rem == 0) = IDS_int(IDS_rem == 0); % Assign integer division where remainder is zero

% Finalize the list of candidate identifiers
CID_FINAL = CIDSS(IDS_final)'; % Extract final CID indices
list.Name1 = repmat(" ", 2000, 1); % Initialize Name1 column
list.Name1(1:length(CID_FINAL)) = SMLS(IDS_final)'; % Fill in with smiles corresponding to CID_FINAL
list.Properties.VariableNames(kSJ2) = FN; % Assign variable names to the list

% Check for duplicates in the final list
[uniqueA, i, j] = unique(CID_FINAL, 'first'); % Find unique CID_FINAL entries
indexToDupes = find(not(ismember(1:numel(CID_FINAL), i))); % Identify duplicates

if ~isempty(indexToDupes) % If duplicates exist
    keyboard % Pause for debugging
end
