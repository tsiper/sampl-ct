images=zeros(761);
for i=1:761
    im_addr=strcat('sheepscan/', num2str(i));
    images(i)= dicomread(im_addr);
end