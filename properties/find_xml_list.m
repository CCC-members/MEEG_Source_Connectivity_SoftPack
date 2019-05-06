function childrens = find_xml_list(file_path,root_tab)
%FIND_XML_PARAMETER Summary of this function goes here
%   Detailed explanation goes here


DOMnode = xmlread(file_path);
element = DOMnode.getElementsByTagName(root_tab);



try
    xml_node = element.item(0);
catch
    result = ["error","The element do not exist."];
    return;
end
childrens = [];
if (xml_node.hasChildNodes)
    childNodes = xml_node.getChildNodes;
    numChildNodes = childNodes.getLength;   
    
    for count = 0: numChildNodes - 1
        child_node = childNodes.item(count);
        if(child_node.getNodeName ~= '#text')
            child = struct();
            child.name = child_node.getNodeName;
            if(child_node.hasAttributes)
                theAttributes = child_node.getAttributes;
                numAttributes = theAttributes.getLength;
                attributes = struct();
                
                for i = 0 : numAttributes - 1
                    attrib = theAttributes.item(i);
                    name = char(attrib.getName);                    
                    value = char(attrib.getValue);
                    attributes.(name)= value;                                
                end
                child.attributes = attributes;
            end
           childrens = [childrens , child];
        end
       
    end
end

end

