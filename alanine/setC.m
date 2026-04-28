function coupling_matrix = setC(coupling_matrix, c1, c2)

    % --- HANDLE 1-BOND COUPLINGS (C1-C2 and C2-C3) ---
    if isscalar(c1)
        j12 = c1;
        j23 = c1;
    else
        j12 = c1(1); % Coupling between Carbon 1 and 2
        j23 = c1(2); % Coupling between Carbon 2 and 3
    end

    % Set C1-C2
    coupling_matrix(1, 2) = j12;
    coupling_matrix(2, 1) = j12;

    % Set C2-C3
    coupling_matrix(2, 3) = j23;
    coupling_matrix(3, 2) = j23;


    % --- HANDLE 2-BOND COUPLINGS (C1-C3) ---
    if isscalar(c2)
        j13 = c2;
    else
        j13 = c2(1);
    end

    % Set C1-C3
    coupling_matrix(1, 3) = j13;
    coupling_matrix(3, 1) = j13;

end