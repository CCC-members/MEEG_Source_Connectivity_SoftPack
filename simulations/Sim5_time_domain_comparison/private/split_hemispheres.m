function [qL,qR,qfull,indvL,indvR,indv,verticesL,verticesR,vertices,facesL,facesR,faces,elec_pos] = split_hemispheres(cortex,elec_pos)
vertices       = cortex.vertices;
center         = min(vertices,[],1) + (1/2)*(max(vertices,[],1) - min(vertices,[],1));
vertices       = vertices - repmat(center,size(vertices,1),1);
elec_pos       = elec_pos - repmat(center,size(elec_pos,1),1);
faces          = cortex.faces;
qfull          = size(vertices,1); 
%% Compute indices and faces correspondign to Left Hemisphere
counter1 = 1; counter2 = 1;indvL = 1;
while counter1
    indfL      = find(sum(ismember(faces,indvL),2));
    indv       = faces(indfL,:);
    indv       = unique([indvL; indv(:)]);
    if isempty(setdiff(indv,indvL));counter1 = 0;
    else indvL = indv;
    end
    counter2   = counter2 + 1;
end
verticesL      = vertices(indvL,:);
facesL         = faces(indfL,:); facesL = facesL - min(facesL(:)) + 1;
qL             = size(verticesL,1);
%%

%% Compute indices and faces correspondign to Right Hemisphere
indvR          = setdiff(1:size(vertices,1),indvL);
indfR          = setdiff(1:size(faces,1),indfL);
verticesR      = vertices(indvR,:);
facesR         = faces(indfR,:); facesR = facesR - min(facesR(:)) + 1;
qR             = size(verticesL,1);
%%