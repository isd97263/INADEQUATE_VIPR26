
alpha_coupling = zeros(6);
alpha_coupling = setC(alpha_coupling, [15 15 15 15 15], [0 0 0 0], [0 0 0]);

beta_coupling = zeros(6);
beta_coupling = setC(beta_coupling, [15 15 15 15 15], [0 0 0 0], [0 0 0]);

s = inad(alpha_coupling, beta_coupling);


