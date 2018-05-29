function [ GridCoordinates ] = ConstructGridCartesian( N, delta, x_offset, y_offset, OverGridFactor)


NewOverGridFactor = ceil(N*OverGridFactor/2)*2/N;
if OverGridFactor ~= NewOverGridFactor
    disp(['Over grid factor was corrected from ',num2str(OverGridFactor),' to ',num2str(NewOverGridFactor)]);
    OverGridFactor = NewOverGridFactor;
end

[ukx,uky] = meshgrid(((1/OverGridFactor:1/OverGridFactor:N)- N/2 -1/OverGridFactor).*delta + x_offset, ((1/OverGridFactor:1/OverGridFactor:N)- N/2 -1/OverGridFactor).*delta + y_offset);

GridCoordinates = [ukx(:) uky(:)];

end

