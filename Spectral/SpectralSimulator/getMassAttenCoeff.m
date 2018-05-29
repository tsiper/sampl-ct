function [ MassAttenCoeff ] = getMassAttenCoeff( material, Energies )
% Calculate the mass attenuation coefficient, for the specific
% material, at the specified energies.
% The material should be from the BodyMaterials struct or a string with 
% chemical formula or simple atomic notation

BodyMaterials = {'adipose', 'blood', 'bone_compact', 'bone_cortical', 'brain', 'lung',...
    'muscle_skeletal', 'muscle_striated', 'skin', 'soft_tissue', 'water'};

if any(cellfun(@(x) isequal(x, material), BodyMaterials))
    MassAttenCoeff = MusFromMaterials(material, Energies);
else
    MassAttenCoeff = XrayMu(material,Energies);
end

end

