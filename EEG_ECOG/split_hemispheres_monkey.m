function [qL,qR,qfull,indvL,indvR,verticesL,verticesR,vertices,facesL,facesR,faces,elec_pos] = split_hemispheres_monkey(cortex,elec_pos)
vertices       = cortex.vertices;
center         = min(vertices,[],1) + (1/2)*(max(vertices,[],1) - min(vertices,[],1));
vertices       = vertices - repmat(center,size(vertices,1),1);
elec_pos       = elec_pos - repmat(center,size(elec_pos,1),1);
faces          = cortex.faces;
qfull          = size(vertices,1); 
%% Compute indices and faces correspondign to Left Hemisphere and Right hemisphere
indvL          = find(vertices(:,1) < 0);
verticesL      = vertices(indvL,:);
qL             = length(indvL);
indvR          = find(vertices(:,1) > 0);
verticesR      = vertices(indvR,:);
qR             = length(indvR);
facesL         = [];
facesR         = [];
%%
cont2L      = 1;
cont2R      = 1;
for cont1 = 1:length(faces)
    if sum(ismember(faces(cont1,:),indvL),2) == 3
        facesL(cont2L,:)    = faces(cont1,:); 
        cont2L              = cont2L + 1;
    elseif sum(ismember(faces(cont1,:),indvR),2) == 3
        facesR(cont2R,:) = faces(cont1,:); 
        cont2R           = cont2R + 1;
    end
end
counter  = 1;
maskL    = ones(size(facesL));
for cont = indvL'
    facesL(facesL == cont*maskL) = counter;
    counter = counter + 1;
end
counter  = 1;
maskR    = ones(size(facesR));
for cont = indvR'
    facesR(facesR == cont*maskR) = counter;
    counter = counter + 1;
end
end
%%