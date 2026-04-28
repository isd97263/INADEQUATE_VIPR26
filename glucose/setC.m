function coupling_matrix = setC(coupling_matrix, c1, c2, c3)
    
    % Get the size of the system from the matrix dimensions
    n = size(coupling_matrix, 1);

    % Set values for 1-bond couplings (J1)
    % Requires a vector c1 of length n-1
    for i = 1:n-1
        coupling_matrix(i, i+1) = c1(i);
        coupling_matrix(i+1, i) = c1(i);
    end

    % Set values for 2-bond couplings (J2)
    % Requires a vector c2 of length n-2
    for i = 1:n-2
        coupling_matrix(i, i+2) = c2(i);
        coupling_matrix(i+2, i) = c2(i);
    end

    % Set values for 3-bond couplings (J3)
    % Requires a vector c3 of length n-3
    for i = 1:n-3
        coupling_matrix(i, i+3) = c3(i);
        coupling_matrix(i+3, i) = c3(i);
    end
end