function testprecondAdjPPFT3
%
% Let D be the preconditioner. Denote by P the 3D ppft and by Q the 3D
% preconditioned ajoint ppft.
%
% This code shows that <A,QB> = <DPA,B>. That is, it shows that we compute
% the preconditioned adjoint ppft correctly.
%
% The printout should be very close to zero.
%
% Yoel Shkolnisky, July 2007.

n=64;
m=3*n+1;
alpha = 2*(n+1)/(n*m);

a=rand(n,n,n);
b=rand(3,3*n+1,n+1,n+1);

ppa=ppft3_ref(a);
ppbT=precondadjppft3_ref(b);

% ppa=ppft3(a);
% ppbT=precondadjppft3(b);


for k=-3*n/2:3*n/2
    for j=-n/2:n/2
        ppa(:,k+3*n/2+1,j+n/2+1,:)=mult(k,alpha,n).*ppa(:,k+3*n/2+1,j+n/2+1,:);
    end
end

for k=-3*n/2:3*n/2
    for j=-n/2:n/2
        ppa(:,k+3*n/2+1,:,j+n/2+1)=mult(k,alpha,n).*ppa(:,k+3*n/2+1,:,j+n/2+1);
    end
end

ip1=ip(ppa,b);
ip2=ip(a,ppbT);

disp(abs(ip1-ip2)/abs(ip1))



% % % % The commented code shows that the function for the preconditioned
% % % % adjoint ppft gives the same answer as applying the preconditioner D
% % % % following by the application of the not-preconditioned adjoint ppft. 
% % % 
% % % n=16;
% % % m=3*n+1;
% % % alpha = 2*(n+1)/(n*m);
% % % 
% % % a=rand(n,n,n);
% % % b=rand(3,3*n+1,n+1,n+1);
% % % ppbT1=ppft3TV3(b);
% % % for k=-3*n/2:3*n/2
% % %     for j=-n/2:n/2
% % %         b(:,k+3*n/2+1,j+n/2+1,:)=mult(k,alpha,n).*b(:,k+3*n/2+1,j+n/2+1,:);
% % %     end
% % % end
% % % 
% % % for k=-3*n/2:3*n/2
% % %     for j=-n/2:n/2
% % %         b(:,k+3*n/2+1,:,j+n/2+1)=mult(k,alpha,n).*b(:,k+3*n/2+1,:,j+n/2+1);
% % %     end
% % % end
% % % ppbT2=adjPPFT3(b);
% % % max(abs(ppbT1(:)-ppbT2(:)))



function v = mult(k,alpha,n)
% Compute preconditioning factor for row k
if k==0
    v = 1/((3*n+1)^2);
else
    v = abs(k*alpha);
end
