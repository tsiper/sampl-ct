function [ ReconstructionParams,I,ChosenMaterials,phantom ] = AddMetal( ReconstructionParams,I,ChosenMaterials,phantom )

numOfMetals = 3;
phantom = [phantom; cell(numOfMetals,1)];
res = ReconstructionParams.PhantomRes;

% Choose metal materials and density:
metal_base = [1 1 1];
MetalMaterials = {'Ti','Fe','Ti'};

% Update scan parameters with metal objects:
ReconstructionParams.base = [ReconstructionParams.base, metal_base];
ChosenMaterials = [ChosenMaterials,MetalMaterials{1:numOfMetals}];

% Choose metal location:
circleCenterX = [50,150, 155]; % to do - adjust to res
circleCenterY = [130, 60, 180]; % to do - adjust to res
circleRadius = [15, 15, 15];

% Add metal to phantom:
x = repmat(1:res,res,1);
y = repmat((1:res)',1,res);
for i=1:numOfMetals
    z = (x-circleCenterX(i)).^2 + (y-circleCenterY(i)).^2;
    phantom{I+i} = zeros(res);
    phantom{I+i}(z <= circleRadius(i)^2) =  ReconstructionParams.base(I+i);
end

I = I+numOfMetals;

end

