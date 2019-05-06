function [] = change_xml_parameter(file_path,root_tab,parameters_name,parameters_value,attributes)
%CHANGE_XML_PARAMETER Summary of this function goes here
%   Detailed explanation goes here

%   file_path           paht of the properties file
%   root_tab            father node of the parameters
%   parameters_name     string vector with all parameter's name
%   parameters_value    string vector with all parameter's value
%   attributes          cell array of key - value of each at
%   cell(['attr1','value1';'attr2','value2'],['attr1','value1'],['attr1','value1';'attr3','value3'])

DOMnode = xmlread(file_path);
bndbox_elem = DOMnode.getElementsByTagName(root_tab);

for i = 1: length(parameters_name)
    element = bndbox_elem.item(0).getElementsByTagName(parameters_name(i));
    if(parameters_value(i) ~= 'null')
        element.item(0).setTextContent(parameters_value(i));
    end
    if(~isempty(attributes))
        attrs_cell = attributes(i);
        attrs_cell = attrs_cell{1,1};
        for m = 1 : 2 :length(attrs_cell)
            attr_name = attrs_cell(m);
            attr_value = attrs_cell(m + 1);
            element.item(0).setAttribute(attr_name,attr_value);
        end
    end
end

xmlwrite(file_path,DOMnode);

end

