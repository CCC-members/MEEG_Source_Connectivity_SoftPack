function result = find_xml_parameter(file_path,root_tab,parameter_name)
%FIND_XML_PARAMETER Summary of this function goes here
%   Detailed explanation goes here


DOMnode = xmlread(file_path);
bndbox_elem = DOMnode.getElementsByTagName(root_tab);
element = bndbox_elem.item(0).getElementsByTagName(parameter_name);

try
    xml_node = element.item(0).getFirstChild;
catch
    result = ["error","The element do not exist."];
    return;
end
if (~isempty(xml_node))
    node_value = char(element.item(0).getFirstChild.getData);
else
    node_value = "null";
end

result = [parameter_name ,  node_value];

attributes = element.item(0).getAttributes;

for i = 0 : attributes.getLength - 1
    attr_name = attributes.item(i).getName;
    attr_value = attributes.item(i).getValue;
    result = [result ;...
        string(attr_name) , string(attr_value)];
end

end

