function [properties,folder] = get_properties(file_path)
% Authors:
% - Deirel Paz Linares
% - Eduardo Gonzalez Moreira
% - Pedro A. Valdes Sosa

% Date: March 16, 2019

% Updates
% - Ariosky Areces Gonzalez

% Date: March 22, 2019

properties = struct;

properties.run_mode = find_xml_parameter(file_path,'properties','run_mode',true);
properties.run_parallel = find_xml_parameter(file_path,'properties','run_parallel',true);
properties.run_frequency_bin = find_xml_parameter(file_path,'properties','run_frequency_bin',true);
properties.run_single_subject = find_xml_parameter(file_path,'properties','run_single_subject',true);
properties.freqres = str2double(find_xml_parameter(file_path,'properties','freq_resol',true));
properties.samplfreq = str2double(find_xml_parameter(file_path,'properties','samp_freq',true));
properties.maxfreq = str2double(find_xml_parameter(file_path,'properties','max_freq',true));


delta_band = find_xml_parameter(file_path,'properties','delta_band',false);
theta_band = find_xml_parameter(file_path,'properties','theta_band',false);
alpha_band = find_xml_parameter(file_path,'properties','alpha_band',false);
beta_band = find_xml_parameter(file_path,'properties','beta_band',false);

frequencies = [];

if(string(delta_band.item(1).getFirstChild.getData) == 'true')
    frequencies = [frequencies;...
        string(delta_band.item(2).getFirstChild.getData),string(delta_band.item(0).getFirstChild.getData), "delta"];
end
if(string(theta_band.item(1).getFirstChild.getData) == 'true')
    frequencies = [frequencies;...
        string(theta_band.item(2).getFirstChild.getData),string(theta_band.item(0).getFirstChild.getData), "theta"];
end
if(string(alpha_band.item(1).getFirstChild.getData) == 'true')
    frequencies = [frequencies;...
        string(alpha_band.item(2).getFirstChild.getData),string(alpha_band.item(0).getFirstChild.getData), "alpha"];
end
if(string(beta_band.item(1).getFirstChild.getData) == 'true')
    frequencies = [frequencies;...
        string(beta_band.item(2).getFirstChild.getData),string(beta_band.item(0).getFirstChild.getData), "beta"];
end

properties.frequencies = frequencies;
folder = find_xml_parameter(file_path,'properties','data_path',true);

end

