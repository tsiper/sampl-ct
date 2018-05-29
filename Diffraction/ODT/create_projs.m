function [ teta_array] = create_projs( RI, lambda, dX, teta, nm )
RI=RI-nm;
[n,~,m]=size(RI);
L1= length(teta);
teta_array= zeros(m, n, L1);

% rotation

%teta- rotation around x

for l=1:L1 %for each angle
    
n_t= zeros(n,n,m);
for t=1:m %For each slice
n_t(t,:, :) = imrotate((squeeze(RI(t,:,:))), teta(l) , 'nearest','crop');
end

teta_array(:,:,l)= 2*pi*dX*sum(n_t,3)/lambda;

end

end

