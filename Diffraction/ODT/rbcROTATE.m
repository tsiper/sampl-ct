function [teta_array] = rbcROTATE( RBC, teta, size, dX, K )

L1= length(teta);
teta_array= zeros(size, size, L1);

% rotation

%teta- rotation around x

for l=1:L1 %for each angle
    
RBC_t= zeros(size,size,size);
for t=1:size %For each slice
RBC_t(t,:, :) = imrotate((squeeze(RBC(t,:,:))), teta(l) , 'nearest','crop');
end

teta_array(:,:,l)= K*dX*sum(RBC_t,3);

end


end


