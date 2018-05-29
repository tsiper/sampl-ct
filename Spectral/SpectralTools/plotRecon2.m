function plotRecon2( x, w, psnr, fval, fig_handel,iter,base)

figure(fig_handel)
subplot(231);imagesc(x{1}); colorbar('EastOutside');
axis off;

subplot(232);imagesc(x{2}); colorbar('EastOutside');
axis off;

subplot(233);loglog(fval(2:end),'LineWidth',1); title(['Function value ',num2str(fval(end))]);
grid on;
subplot(236);semilogx(psnr(2:end),'LineWidth',1);
title(['PSNR value ',num2str(psnr(end))]);
grid on;

subplot(234);imagesc(w{1}/base(1)); colorbar('EastOutside');
subplot(235);imagesc(w{2}/base(2)); colorbar('EastOutside');

suptitle(['FISTA, Iteration ',num2str(iter)]);
drawnow;

end

