function [ x_star ] = FGP( b, lambda , N_Iters )
%FGP Algorithm according to Amir Beck's
% Solves the following denoising problem
% x_star = arg min_b' ||b'-b||_2^2+lambda*TV(b')
% The algorith runs for N iterations - usually 20-30 is more than enough


% Getting the dimensions of the objective
[m,n] = size(b);

% Initializing the step size to 1
t_k = 1;

% Generating the initial derivatives, r1 and s1
r_k = zeros(m-1,n); p_k_prev = r_k;
s_k = zeros(m,n-1); q_k_prev = s_k;

% Running the main iteration of the FGP algorithm
for k=1:N_Iters
    
    % Getting the inner arguments for the projection (eq 4.9)
    [p_temp,q_temp] = Lt( Proj_C( b - lambda*L(r_k,s_k) ) );

    p_temp = r_k + 1/(8*lambda) * p_temp;
    q_temp = s_k + 1/(8*lambda) * q_temp;
    
    [p_k , q_k] = Proj_P(p_temp , q_temp);
    
    % Updating step size
    t_k_next = (1+ sqrt(1+4*t_k^2)) / 2;
    
    % Computing next derivatives
    r_k = p_k + (t_k-1)/t_k_next * (p_k - p_k_prev);
    s_k = q_k + (t_k-1)/t_k_next * (q_k - q_k_prev);
    
    % Setting variables towards next iteration
    t_k = t_k_next;
    p_k_prev = p_k;
    q_k_prev = q_k;
    
end

x_star = Proj_C(b - lambda * L (p_k,q_k));

end