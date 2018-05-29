%Author: Stefano "grog" Ghio, Michele Bozzano, Bruno Mazzarello
%Released under CreativeCommons and GPLv2 licenses

%EM - Expectation Maximization
%Formula is:
%f_k+1 = (f_k / alpha) (At (g / (A f_k)))
%f_k is solution at k-th iteration, at first iteration it is our guess
%g is image sinogram
%A f_k is Radon transform of f_k
%At is inverse Radon of its argument
%alpha is inverse Radon of a sinogram with all values 1 which represents our scanning machine

%clear matlab environment
clear; close all;
N = 256;
theta=0:135; %limited projection angle 90 degrees instead of 180
F=phantom(N); %create test phantom 64x64 pixels aka alien
figure, imshow(F,[]),title('Alien vs Predator'); %show alien
S1 = sum(sum(F));%calculate pixels sum on F
R = radon(F,theta); %apply Radon aka scan the alien to create a sinogram
% figure, imagesc(R),title('Sinoalien');%show sinogram
%set all values <0 to 0
index=R<0;
R(index)=0;
n = 2000;%iterations
Fk=ones(N);%our initial guess
%create alpha aka pixels sum of the projection matrix (our scanning machine)
% sinogramma1=ones(N,length(theta));
sinogramma1=ones(round(sqrt(2)*(N+1))+3,length(theta));
alpha=iradon(sinogramma1, theta,'linear', 'none', 1,N);
%calculate relative error
relerrfact=sum(sum(F.^2));
wb = MyWaitbar(0,'EM Reconstruction');
for  k=1:n
     Afk = radon(Fk,theta);%create sinogram from current step solution
     %calculate g / A f_k=Noised./(A f_k+eps); aka initial sinogram / current step sinogram
     %eps is matlab thing to prevent 0 division
     GsuAFK=R./(Afk+eps);
     retro=iradon(GsuAFK, theta, 'linear', 'none', 1,N);%At (g / (A f_k))
     %multiply current step with previous step result and divide for alpha updating f_k
     ratio=Fk./alpha;
     Fk=ratio.*retro;
     %normalize
     St = sum(sum(Fk));
     Fk = (Fk/St)*S1;
	 %calculate step improvement
     Arrerr(k) = sum(sum((F - Fk).^2))/relerrfact;
	 %stop when there is no more improvement
     if((k>2) &&(Arrerr(k)>Arrerr(k-1)))
        break;
     end
     MyWaitbar(k/n);
end
close(wb);
% figure, imagesc(Fk),title('Fk');%show reconstructed alien
% figure, plot(Arrerr),title('Arrerr');%show error improvement over all iterations

%compare our result with the one we would have had using the FBP - Filtered Back Projection
% easy=iradon(R,theta, 'Shepp-Logan',1,N);
% figure, imagesc(easy),title('FBP');

%calculate error between EM and FBP - with limited image size and projection degree FBP is bad!
% FBPerr=sum(sum((F - easy).^2))/relerrfact;

%% Measuring results

PSNR_EM  = psnr(Fk,F);
PSNR_FBP = psnr(easy,F);

%% Plotting

figure;
subplot(221);imagesc(F);title('Ground Truth');
subplot(222);imagesc(easy);title(['FBP Reconstruction ',num2str(PSNR_FBP)]);
subplot(223);imagesc(Fk);title(['EM Reconstruction ',num2str(PSNR_EM)]);
subplot(224);plot(Arrerr);title('Solution Convergence');
