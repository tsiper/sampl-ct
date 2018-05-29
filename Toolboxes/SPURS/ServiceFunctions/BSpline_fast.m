function x=BSpline_fast(t,SplineDegree)

[tm, tn] = size(t);
t = t(:);
n = SplineDegree;
rel_t_idx = find(abs(t)<SplineDegree);

if isempty(rel_t_idx)
    x = zeros(tm,tn);
else
    sdegree = (0:(n+1));
    NCKv = zeros(n+2,1);
    for ii=0:n+1
        NCKv(ii+1) = nchoosek(n+1,ii);
    end
    NCK_1_K = repmat(NCKv.*(-1).^(sdegree.'),1,length(rel_t_idx));
    X_K_n_1_2 = (repmat(t(rel_t_idx).',length(NCKv),1)-repmat(sdegree.',1,length(rel_t_idx))+(n+1)/2);
    X_K_n_1_2 = (X_K_n_1_2 > 0).*X_K_n_1_2;
    X_K_n_1_2 = X_K_n_1_2.^n;
    B = NCK_1_K.*X_K_n_1_2;
    xt = zeros(size(t));
    xt(rel_t_idx)=(1/factorial(n)).*sum(B,1);
    x = reshape(xt,tm,tn);
end
end

