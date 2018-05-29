function [x_min] = conjgrad(A,b,x,iters,eps)
% x = x(:);

if nargin < 5
    eps = 1e-8;
end
if nargin < 4
    iters = 20;
end

% r=b-A(vec2im(x));
r=b-A(x);
p=r;
rs_old = r(:)'*r(:);
rs_min = rs_old;
x_min  = zeros(size(x));
for i=1:iters
    Ap = A(p);

    alpha=rs_old/(p(:)'*Ap(:));
    x=x+alpha*p;
    r=r-alpha*Ap;
    rs_new=r(:)'*r(:);


    if (sqrt(rs_new)<eps)
%         disp('ConjGrad Stopped - Error threshold met');
        break;
    end
    p=r+(rs_new/rs_old)*p;
    rs_old=rs_new;
%     display(sprintf('Conj-Grad, Iteration %d, Err %g\n',i,norm(rs_new)));
%     if PlotFlag
%         subplot(224);
%         imagesc(A(vec2im(x)));colorbar;
%         title('Sinogram');
%         subplot(223);
%         imagesc(vec2im(x)); ssim
%         title('Solution');
%         drawnow;
%     end
    if rs_new<rs_min
        x_min = x;
        rs_min = rs_new;
    end
end
end