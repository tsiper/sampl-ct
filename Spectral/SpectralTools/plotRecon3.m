function plotRecon3( x, w, psnr, fval, fig_handel,iter,base)

figure(fig_handel)
subplot(231);imagesc(x{1}); colorbar('EastOutside');
axis off;
subplot(232);imagesc(x{2}); colorbar('EastOutside');
axis off;
subplot(233);imagesc(x{3}); colorbar('EastOutside');
axis off;
subplot(234);loglog(fval(2:end),'LineWidth',1); title(['Function value ',num2str(fval(end))]);
grid on;
subplot(235);semilogx(psnr(2:end),'LineWidth',1);
title(['PSNR value ',num2str(psnr(end))]);
grid on;

alpha_hat(:,:,1) = w{1}/base(1);
alpha_hat(:,:,2) = w{2}/base(2);
alpha_hat(:,:,3) = w{3}/base(3);
alpha_hat = (alpha_hat-min(alpha_hat(:)))/range(alpha_hat(:));
subplot(236);imagesc(alpha_hat);
suptitle(['FISTA, Iteration ',num2str(iter)]);
drawnow;

end

