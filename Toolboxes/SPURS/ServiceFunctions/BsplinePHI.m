function PHI=BsplinePHI(ku, knu,spline_degree);
% function A=interp2spline(xyc, xyk)
% 1D interpolation using spline function.
% ku: the column vector of Cartesian k-space positions,
% knu: the column vector of Non-Cartesian k-space postion.
% PHI: the linear transform matrix s.t. the (DATA@knu)=PHI*(DATA@ku);


length_ku = length(ku);
length_knu = length(knu);

disp(['Size of PHI matrix is (',num2str(length_knu),':',num2str(length_ku),')']);

max_size = 3000;
% A = spalloc(length_xnu,length_xu,round(length_xnu*length_xnu*0.001));
ExpectedNNZ = round(length_knu*pi()*((1/abs(ku(1)-ku(2)))*spline_degree)^2);
disp(['PHI is expected to have at most ',num2str(ExpectedNNZ),' non-zeros']);
disp('Allocating memory for PHI');

i = zeros(ExpectedNNZ,1);
j = zeros(ExpectedNNZ,1);
s = zeros(ExpectedNNZ,1);

iimax = ceil(length_knu/max_size);
jjmax = ceil(length_ku/max_size);
h = waitbar(0,'Calculating PHI matrix...');
lv = 0;
for ii=1:iimax
    iim = (ii-1)*max_size+1;
    iin = min(ii*max_size,length_knu);
    rxnu = knu(iim:iin);
    for jj=1:jjmax
        jjm = (jj-1)*max_size+1;
        jjn = min(jj*max_size,length_ku);
        rxu = ku(jjm:jjn);
        pxu_m_pxnu = repmat(rxu.',length(rxnu),1)-repmat(rxnu,1,length(rxu));
        Atemp=BSpline_fast(real(pxu_m_pxnu),spline_degree).*BSpline_fast(imag(pxu_m_pxnu),spline_degree);
        [itemp,jtemp,stemp] = find(Atemp);
        if ~isempty(itemp)
            i(lv+1:lv+length(itemp)) = itemp + iim-1;
            j(lv+1:lv+length(jtemp)) = jtemp + jjm-1;
            s(lv+1:lv+length(stemp)) = stemp;
            lv = lv + length(itemp);
        end
        waitbar(((ii-1)*jjmax+jj) / (iimax*jjmax))
    end
end
close(h)
PHI = sparse(i(1:lv),j(1:lv),s(1:lv),length_knu,length_ku);
disp(['PHI matrix has ',num2str(nnz(PHI)),' non-zeros']);
disp(['PHI matrix is ',num2str(100-nnz(PHI)/length_ku/length_knu*100),'% sparse']);

