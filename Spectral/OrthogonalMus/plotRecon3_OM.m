function plotRecon3_OM( x, w, psnr, fval, fig_handel,iter,invP_operator,base )

figure(fig_handel)
subplot(331);imagesc(x{1}); colorbar('EastOutside');
axis off;
subplot(332);imagesc(x{2}); colorbar('EastOutside');
axis off;
subplot(333);imagesc(x{3}); colorbar('EastOutside');
axis off;
PhantomRes = size(w{1},1);
RadonRes = size(x{1}); 
I = length(w);
column_stack_w = zeros(PhantomRes*PhantomRes,I);
column_stack_a = zeros(RadonRes(1)*RadonRes(2),I);
for m=1:I
    column_stack_w(:,m) = w{m}(:);
    column_stack_a(:,m) = x{m}(:);
end
column_stack_true_w = (invP_operator*column_stack_w')';
column_stack_true_a = (invP_operator*column_stack_a')';
alpha_hat= zeros(PhantomRes,PhantomRes,I);
a_hat = zeros([RadonRes,I]);
for m=1:I
    alpha_hat(:,:,m) = reshape(column_stack_true_w(:,m),[PhantomRes,PhantomRes])/base(m);
    a_hat(:,:,m) = reshape(column_stack_true_a(:,m),RadonRes);
end
alpha_hat = (alpha_hat-min(alpha_hat(:)))/range(alpha_hat(:));

subplot(334);imagesc(a_hat(:,:,1)); colorbar('EastOutside');
axis off;
subplot(335);imagesc(a_hat(:,:,2)); colorbar('EastOutside');
axis off;
subplot(336);imagesc(a_hat(:,:,3)); colorbar('EastOutside');
axis off;

subplot(337);loglog(fval(2:end),'LineWidth',1);
title(['Function value ',num2str(fval(end))]);
grid on;
subplot(338);semilogx(psnr(2:end),'LineWidth',1);
title(['PSNR value ',num2str(psnr(end))]);
grid on;

subplot(339);imagesc(alpha_hat);
suptitle(['FISTA, Iteration ',num2str(iter)]);
drawnow;

end

