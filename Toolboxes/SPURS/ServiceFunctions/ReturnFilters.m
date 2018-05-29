function [ r_aq,ImgFilter ] = ReturnFilters( KernelFunction,OverGridFactor,sqrtN )

r_aq = BSpline_fast((-ceil(KernelFunction.degree/2*OverGridFactor):1:ceil(KernelFunction.degree/2*OverGridFactor))./OverGridFactor,KernelFunction.degree);
r_aq = (r_aq + fliplr(r_aq))./2;
r_aq = r_aq./sum(r_aq);

xx = (((-sqrtN/2:1/OverGridFactor:(sqrtN/2-1/OverGridFactor))./sqrtN)*OverGridFactor);
R_AQ = (sinc(xx)).^(KernelFunction.degree+1);
R_AA = (xx>-1/2).*(xx<1/2);
ImgFilter = repmat(R_AQ.*R_AA,sqrtN*OverGridFactor,1);
ImgFilter = ImgFilter.*(ImgFilter');

end

