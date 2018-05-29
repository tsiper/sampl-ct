function [X,info] = symkaczmarz(A,b,K,x0,options)
%SYMKACZMARZ Symmetric Kaczmarz method
%
%   [X,info] = symkaczmarz(A,b,K)
%   [X,info] = symkaczmarz(A,b,K,x0)
%   [X,info] = symkaczmarz(A,b,K,x0,options)
%
% Implements the symmetric Kaczmarz method: each sweep consists of a Kaczmarz
% sweep followed by a Kaczmarz sweep with the equations in reverse order.
%
% Input:
%   A        m times n matrix.
%   b        m times 1 vector.
%   K        Number of iteration. If K is a scalar, then K is the maximum
%            number of iterations and only the last iterate is saved.
%            If K is a vector, then the largest value in K is the maximum
%            number of iterations and only iterates corresponding to the
%            values in K are saved, together with the last iterate.
%            If K is empty then a stopping criterion must be specified.
%   x0       n times 1 starting vector. Default: x0 = 0.
%   options  Struct with the following fields:
%       lambda    The relaxation parameter. If lambda is a scalar then
%                 the corresponding value is used in each iteration; the
%                 default value is 1.
%                 If lambda is a string, then it refers to a method to 
%                 determine lambda in each iteration. For this method the
%                 following strings can be specified:
%                     'psi1' : lambda is chosen using the Psi_1-based 
%                                 relaxation method.
%                     'psi2' : lambda is chosen using the Psi_2-based
%                                 relaxation method.
%       stoprule  Struct containing the following information about the
%                 stopping rule:
%                     type = 'none' : (Default) the only stopping rule
%                                     is the maximum number of iterations.
%                            'NCP'  : Normalized Cumulative Peridogram.
%                            'DP'   : Discrepancy Principle.
%                     taudelta = the product of tau and delta, only
%                                necessary for DP.
%       nonneg    Logical; if true then nonnegativity in enforced in
%                 each step.
%       box       Upper bound L in box constraint [0,L] on pixel values.
%       damping   A parameter D to avoid division by very small row norms
%                 by adding D*max_i{||a^i||_2^2} to ||a^i||_2^2.
%
% Output:
%   X     Matrix containing the saved iterations.
%   info  Information vector with 2 elements.
%         info(1) = 0 : stopped by maximum number of iterations
%                   1 : stopped by NCP-rule
%                   2 : stopped by DP-rule
%         info(2) = no. of iterations.
%
% See also: kaczmarz, randkaczmarz

% Maria Saxild-Hansen and Per Chr. Hansen, July 5, 2015, DTU Compute.

% Reference: Å. Björck and T. Elfving, Accelerated projection methods for 
% computiong pseudoinverse solutions of systems of linear equations, BIT 
% 19 (1979), pp. 145-163.

[m,n] = size(A);
A = A';  % Faster to perform sparse column operations.

% Check the number of inputs.
if nargin < 3
    error('Too few input arguments')
end

% Default value of starting vector x0.
if nargin < 4 || isempty(x0)
    x0 = zeros(n,1);
end

% Input check: The sizes of A, b and x must match.
if size(b,1) ~= m || size(b,2) ~= 1
    error('The size of A and b do not match')
elseif size(x0,1) ~= n || size(x0,2) ~= 1
    error('The size of x0 does not match the problem')
end

% Initialization.
if nargin < 5
    if isempty(K)
        error('No stopping rule specified')
    end
    stoprule = 'NO';
    lambda = 1;
    Knew = sort(K);
    kmax = Knew(end);
    X = zeros(n,length(K));
    casel = 1;
    
    % Default is no nonnegativity or box constraint or damping.
    nonneg = false;
    boxcon = false;
    damp = 0;

else
% Check the contents of options, if present.
    
    % Nonnegativity.
    if isfield(options,'nonneg')
        nonneg = options.nonneg;
    else
        nonneg = false;
    end
    
    % Box constraints [0,L].
    if isfield(options,'box')
        nonneg = true;
        boxcon = true;
        L = options.box;
    else
        boxcon = false;
    end
    
    % Damping.
    if isfield(options,'damping')
        damp = options.damping;
        if damp<0, error('Damping must be positive'), end
    else
        damp = 0;
    end
   
    if isfield(options,'stoprule') && isfield(options.stoprule,'type')
        stoprule = options.stoprule.type;
        if ischar(stoprule)
            if strncmpi(stoprule,'DP',2)
                % DP stopping rule.
                if isfield(options.stoprule,'taudelta')
                    taudelta = options.stoprule.taudelta;
                else
                    error('The factor taudelta must be specified when using DP')
                end   
                
                % Check that the first iteration should be performed:
                rk = b - A'*x0;  % Remember that A is transposed.
                nrk = norm(rk);
                
                if nrk <= taudelta
                    info = [2 0];
                    X = x0;
                    return
                end % end the DP-rule.
                
            elseif strncmpi(stoprule,'NC',2)
                % NCP stopping rule.
                dk = inf;
                q = floor(m/2);
                c_white = (1:q)'./q;
                
                if ~isempty(K)
                    K = [K max(K)+1];
                end
            elseif strncmpi(stoprule,'NO',2)
                % No stopping rule.
                if isempty(K)
                    error('No stopping rule specified')
                end
                
            else
                % Other stopping rules.
                error('The shosen stopping rule is not valid')
            end % end different stopping rules
        else
            error('The stoprule type must be a string')
        end % end stoprule is a string.
        if isempty(K)
            kmax = inf;
            X = zeros(n,1);
        else
            Knew = sort(K);
            kmax = Knew(end);
            X = zeros(n,length(K));
        end
    else
        if isempty(K)
            error('No stopping rule specified')
        else
            Knew = sort(K);
            kmax = Knew(end);
            X = zeros(n,length(K));
            stoprule = 'NO';
        end
    end % end stoprule type specified.
    
    if isfield(options,'lambda')
        lambda = options.lambda;
        % If lambda is a scalar.
        if ~ischar(lambda)
            
            % Convergence check.
            if lambda < 0 || lambda > 2
                warning('MATLAB:UnstableRelaxParam',...
                    'The lambda value is outside the interval (0,2)');
            end
            casel = 1;
        else
            % Calculates the lambda value according to the chosen method.
            if strncmpi(lambda,'psi1',4)
                % Method: Psi1
                casel = 3;
                ite = 0;
                
                sigma1tildesquare = 1;
                       
                % Precalculates the roots for the case where the stopping
                % iteration is unknown.
                if isempty(K)
                    % Precalculate z for 1000 iterations.
                    z = calczeta(2:999);
                else
                    % Precalculate the roots for the known number of
                    % iterations.
                    z = calczeta(2:max(K)-1);
                end
                
                % Define the values for lambda according to Psi1.

                lambdak = [sqrt(2); sqrt(2); 2*(1-z)]/sigma1tildesquare;
                
            elseif strncmpi(lambda,'psi2',4)
                % Method: Psi2.
                casel = 3;
                ite = 0;
                
                sigma1tildesquare = 1;
                
                % Precalculates the roots for the case where the stopping
                % iteration is unknown.
                if isempty(K)
                    % Precalculate z for 1000 iterations.
                    kk = 2:999;
                    z = calczeta(kk);
                else
                    % Precalculates the roots for the known number of
                    % iteerations.
                    kk = 2:max(K)-1;
                    z = calczeta(kk);
                end
                
                % Define the values for lambda according to the psi2.
                lambdak = [sqrt(2); sqrt(2); 
                    2*(1-z)./((1-z.^(kk')).^2)]/sigma1tildesquare;
            else
                error(['The chosen relaxation strategy is not ',...
                    'valid for this metod.'])
            end % end check of lambda strategies.
        end
        
    else
        casel = 1;
        lambda = 1;
    end
end % end if nargin includes options.

% Initialization before iterations.
xk = x0;
normAi = full(abs(sum(A.*A,1)));  % Remember that A is transposed.
I = find(normAi>0);
I = [I, I(end-1:-1:2)];
normAi = normAi + damp*max(normAi);

stop = 0;
k = 0;
l = 0;
klast = 0;

while ~stop
    k = k + 1;
    if strncmpi(stoprule,'NC',2)
        xkm1 = xk;
    end
    % The Kaczmarz sweep followed by the reverse order sweep.
    if casel == 1
        for i = I
            ai = full(A(:,i));  % Remember that A is transposed.
            xk = xk + (lambda*(b(i) - ai'*xk)/normAi(i))*ai;
            if nonneg, xk = max(xk,0); end
            if boxcon, xk = min(xk,L); end
        end

    elseif casel == 3
        ite = ite + 1;
        % Check if we need to precalculate another 1000 values of lambda.
        if isempty(K)
            if k ~= 1 && rem(k-1,1000) == 0
                ite = 1;
                kk = k:k+999;
                z = calczeta(kk);
                if strncmpi(lambda,'psi1',4)
                    lambdak = (2*(1-z))/sigma1tildesquare;
                else
                    lambdak = (2*(1-z)./((1-z.^(kk')).^2))/sigma1tildesquare;
                end
            end
        end
        for i = I
            ai = full(A(:,i))';  % Remember that A is transposed.
            xk = xk + (lambdak(ite)*(b(i) - ai*xk)/normAi(i))*ai';
        end
    end

    % Stopping rules:
    if strncmpi(stoprule,'DP',2)
        % DP stopping rule.
        nrk = norm(b-A'*xk);  % Remember that A is transposed.
        if nrk <= taudelta || k >= kmax
            stop = 1;
            if k ~= kmax
                info = [2 k];
            else
                info = [0 k];
            end
        end
        
    elseif strncmpi(stoprule,'NC',2)
        % NCP stopping rule.
        rkh = fft( b-A'*xk);  % Remember that A is transposed.
        pk = abs(rkh(1:q+1)).^2;
        c = zeros(q,1);
        for index = 1:q
            c(index) = sum(pk(2:index+1))/sum(pk(2:end));
        end
        if dk < norm(c-c_white) || k >= kmax
            stop = 1;
            if k ~= kmax
                info = [1 k-1];
            else
                info = [0 k-1];
            end
        else
            dk = norm(c-c_white);
        end
    elseif strncmpi(stoprule,'NO',2)
        % No stopping rule.
        if k >= kmax
            stop = 1;
            info = [0 k];
        end
    end % end stoprule type.
    
    % If the current iteration is requested saved.
    if (~isempty(K) && k == Knew(l+1)) || stop
        l = l + 1;
        % Savs the current iteration.
        if strncmpi(stoprule,'NC',2)
            if ~(stop && klast == k-1)
                X(:,l) = xkm1;
            else
                l = l - 1;
            end
        else
            X(:,l) = xk;
        end
        klast = k;
    end
    
end
X = X(:,1:l);