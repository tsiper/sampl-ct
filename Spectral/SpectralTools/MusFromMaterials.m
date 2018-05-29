function mus = MusFromMaterials(materials,egys)
%    function mus = MusFromMaterials(materials,egys)
%    make the mus for materials -- limit to body materials now
%    14/May/2014 9:59
%    COPYRIGHT: © Robert E. Alvarez, 2014. All rights reserved.

if iscell(materials)
    for k = 1:length(materials)
        assert(ischar(materials{k}));
    end
elseif ischar(materials)
    materials = {materials}; % make it into cell array
else
    error('materials has to be string or cell array of strings');
end
    
  % test for 1d VECTOR
assert(isa(egys,'double') && isvector(egys));
   % the mus for 'soft_tissue','bone_cortical'
bdat = BodyMaterialCompositionFunc;  % load the composition tables
%  fields of bdat are: 'adipose','blood','bone_compact','bone_cortical','brain','lung','muscle_skeletal','muscle_striated','skin','soft_tissue' ...
assert(all(isfield(bdat,materials)));

nmaterials = length(materials);
mus = zeros(numel(egys),nmaterials);
for kmaterial = 1:nmaterials
  dz = bdat.(char(materials{kmaterial})); % bdat is struct from BodyMaterialCompositionFunc so this a ztable
  mus(:,kmaterial) = XrayMu(dz(2:end,:),egys)'; % the first row of dz has the density and is not used by xraymu
end
