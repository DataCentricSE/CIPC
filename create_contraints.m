function [A,b,Aeq,beq] = create_contraints(CID,u)
    % Reshape CID into a column vector to ensure proper dimensions for processing
    CID = reshape(CID',[],1)'; 

    % Identify unique values in CID while preserving their order
    [~, w] = unique(CID, 'stable'); 
    duplicate_indices = setdiff(1:numel(CID), w); % Find indices of duplicate entries
    duplicate_values = unique(CID(1, duplicate_indices), 'stable'); % Extract unique duplicate values

    % Initialize the inequality constraint matrix A
    A  = zeros(numel(duplicate_values), numel(CID).^2);
    
    % Loop through each unique duplicate value to fill in the constraint matrix A
    for k_A = 1:size(duplicate_values,2)
        ID_ones = duplicate_values(k_A) == CID; % Find indices where current duplicate value matches CID
        nums = find(ID_ones); % Get the positions of these indices
        A_temp = zeros(numel(CID), numel(CID)); % Temporary matrix for the current duplicate
        
        % Populate the temporary matrix to mark that only one of the candidates can be selected
        for k_A2 = 1:sum(ID_ones)
            A_temp(nums(k_A2), nums(k_A2)) = 1; % Set diagonal element for this candidate
        end
        A(k_A,:) = reshape(A_temp', 1, []); % Reshape the temporary matrix into a row and store in A
    end

    % Initialize the equality constraint matrix A_eq2 for each unique peak
    u_u = unique(u); 
    A_eq2 = zeros(numel(u_u), numel(u)^2);

    % Loop through each unique peak to create the equality constraints
    for k = 1:length(unique(u))
        A_eq_temp = zeros(numel(u), numel(u)); % Temporary matrix for the current peak
        A_eq_part = repmat(u == u_u(k), sum(u == u_u(k)), 1); % Identify candidates corresponding to this peak
        indexes = find(u == u_u(k)); % Get indices of these candidates
        A_eq_temp(indexes, :) = A_eq_part; % Fill the temporary matrix with matching candidates
        A_eq2(k,:) = reshape(A_eq_temp', 1, []); % Reshape and store the equality matrix
    end

    % Initialize the equality constraint matrix A_eq
    A_eq = zeros(numel(u), numel(u)^2);

    % Populate the equality constraint matrix A_eq
    for k = 1:length(u)
        A_eq_temp = zeros(numel(u), numel(u)); % Create a temporary matrix
        A_eq_temp(k,:) = -1; % Set the row to -1 for the current candidate
        A_eq_temp(:, k) = 1; % Set the column to 1 for the current candidate
        A_eq_temp(k,k) = 0; % Set the diagonal element to 0 (to prevent counting itself)
        A_eq(k,:) = reshape(A_eq_temp', 1, []); % Reshape and store in A_eq
    end

    % Combine the equality constraints
    Aeq = [A_eq; A_eq2; reshape(eye(numel(u))', 1, []); ones(1, numel(u)^2)];

    % Right-hand side for equality constraints
    beq = [zeros(numel(u), 1); ones(numel(u_u), 1); numel(u_u); numel(u)]; 

    % Right-hand side for inequality constraints, ensuring candidates do not exceed selection
    b = ones(size(A, 1), 1); 
end
