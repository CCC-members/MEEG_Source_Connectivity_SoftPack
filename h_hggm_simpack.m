classdef h_hggm_simpack < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        HHGGMSimulationPackUIFigure   matlab.ui.Figure
        FileMenu                      matlab.ui.container.Menu
        SelectdatafileMenu            matlab.ui.container.Menu
        SelectoutputfolderMenu        matlab.ui.container.Menu
        ExitMenu                      matlab.ui.container.Menu
        ToolsMenu                     matlab.ui.container.Menu
        HGGMMenu                      matlab.ui.container.Menu
        HGGMConvergenceJankovaConditionsMenu  matlab.ui.container.Menu
        HHGGMSimplifiedHeadModelMenu  matlab.ui.container.Menu
        HHGGMRealisticHeadModelMenu   matlab.ui.container.Menu
        HHGGMHeadModelComparisonMenu  matlab.ui.container.Menu
        SSVEPMenu                     matlab.ui.container.Menu
        EEGECOGMenu                   matlab.ui.container.Menu
        AllAnalyzeMenu                matlab.ui.container.Menu
        HelpMenu                      matlab.ui.container.Menu
    end

    
    properties (Access = public)
        output_sourse % Description
        test_data_sourse
        downloaded
        data_url
        count
    end
    
    methods (Access = private)
        
        function result = get_test_data(app,simuling)
            app.data_url = [
                "https://drive.google.com/uc?id=16mBEj8ga90mBvZnWTCtDT-LzH_041E-6",...
                "https://drive.google.com/uc?id=1BdNvxeCOYSmB5jZt1_ueWcbAaWkSs9HP",...
                "https://drive.google.com/uc?id=10gCGzTQj0LVVjUk29TmqL5QP_Dit21BY",...
                "https://drive.google.com/uc?id=1WcFJ07dLDPRAVK4vxWIuMBH8dw5dXnGr",...
                "https://drive.google.com/uc?id=1onnwMmpqLPF5VH_qYi3PHtSjvim6tlYH",...
                "https://drive.google.com/uc?id=1qDIifcmFvdI8bTEQCL89cVGIh3qpP84d",...
                "https://drive.google.com/uc?id=1pk_6KbwEQBqxIq2WMW6ekB_raK3cKrE_",...
                "https://drive.google.com/uc?id=1yPhC_0ueGJO-Xk55ui1RZ3sNah30Nm2c",...
                "https://drive.google.com/uc?id=1uQv9cKB4tjfutxKDkgCfcoZUZb-fMQQN",...
                "https://drive.google.com/uc?id=13R6zqVOHs330gYV8uhjUAunyVBKPXHwR",...
                "https://drive.google.com/uc?id=1DwJdK3Iulu0fv6urMuCjXmSrvAaQv-0q",...
                "https://drive.google.com/uc?id=14sQfnBUUo1bVOhujT8Hn3dtMwiUlvmVi",...
                "https://drive.google.com/uc?id=148wgbdiENxC5hurgiR3RHZsQm411OKrn",...
                "https://drive.google.com/uc?id=1vNCgYYjY-yzmFYMuZUbTlESUUBOd6-uk",...
                "https://drive.google.com/uc?id=11yXPfWNwgC9pvUqbRRxJBQ8bFMiyda94",...
                "https://drive.google.com/uc?id=1eQtQqjUTP5yng8zfXyR6Smg6wpSUz18e"
                ];
            
            result = false;
            file_path = strcat('properties',filesep,'properties.xml');
            root_tab =  'properties';
            parameter_name = "data_downloaded";
            app.count = find_xml_parameter(file_path,root_tab,parameter_name);
            app.count = string(app.count(1,2));
            if(app.count ~= "end")
                answer = questdlg('Did you download the data?', ...
                    'Select data', ...
                    'Yes I did','Download','Cancel','Close');
                % Handle response
                switch answer
                    case 'Yes I did'
                        [file,path] = uigetfile('*.zip');
                        if isequal(file,0)
                            disp('User selected Cancel');
                            return;
                        else
                            disp(['User selected ', fullfile(path,file)]);
                            f = dialog('Position',[300 300 250 80]);
                            
                            iconsClassName = 'com.mathworks.widgets.BusyAffordance$AffordanceSize';
                            iconsSizeEnums = javaMethod('values',iconsClassName);
                            SIZE_32x32 = iconsSizeEnums(2);  % (1) = 16x16,  (2) = 32x32
                            jObj = com.mathworks.widgets.BusyAffordance(SIZE_32x32, 'Unpacking test data...');  % icon, label
                            
                            jObj.setPaintsWhenStopped(true);  % default = false
                            jObj.useWhiteDots(false);         % default = false (true is good for dark backgrounds)
                            javacomponent(jObj.getComponent, [50,10,150,80], f);
                            jObj.start;
                            pause(1);
                            try
                                exampleFiles = unzip(fullfile(path,file),pwd);
                            catch
                                delete(f);
                                errordlg('Unpackage error!!!','Error');
                                return;
                            end
                            jObj.stop;
                            jObj.setBusyText('All done!');
                            pause(2);
                            delete(f);
                            msgbox('Completed unpackage!!!','Info');
                            first_file = strcat('EEG_ECOG',filesep,'data',filesep,'4-Head_Model',filesep,'Electrodes',filesep,'EcoG-elecs_Su.mat');
                            if(isfile(first_file))
                                change_xml_parameter(file_path,root_tab,[parameter_name],["end"],cell(0,0));
                                result = true;
                            else
                                errordlg('The selected package is not correct!!','Error');
                            end
                        end
                    case 'Download'
                        f = dialog('Position',[300 300 250 80]);
                        
                        iconsClassName = 'com.mathworks.widgets.BusyAffordance$AffordanceSize';
                        iconsSizeEnums = javaMethod('values',iconsClassName);
                        SIZE_32x32 = iconsSizeEnums(2);  % (1) = 16x16,  (2) = 32x32
                        jObj = com.mathworks.widgets.BusyAffordance(SIZE_32x32, 'Downloading test data...');  % icon, label
                        
                        jObj.setPaintsWhenStopped(true);  % default = false
                        jObj.useWhiteDots(false);         % default = false (true is good for dark backgrounds)
                        javacomponent(jObj.getComponent, [50,10,150,80], f);
                        jObj.start;
                        pause(1);
                        try
                            app.count = str2double(app.count);
                            for i= app.count : length(app.data_url)
                                disp(strcat("Downloading file ",string(i),"  of ", string(length(app.data_url)) ));
                                jObj.setBusyText(strcat("Downloading file ",string(i),"  of ", string(length(app.data_url)) ));
                                url =  app.data_url(i);
                                last = string(i);
                                if(i<10)
                                    last = strcat('0',string(i));
                                end
                                filename = strcat('H_HGGM_test_data.z',last);
                                if(i == length(app.data_url))
                                    filename = strcat('H_HGGM_test_data.zip');
                                end
                                options = weboptions('Timeout',Inf,'RequestMethod','get');
                                
                                % ------Downloding the .zip file  ----
                                outfilename = websave(filename,url,options);
                                s = dir(char(filename));
                                pause(1);
                                
                                % ----- Checking if the .zip file is correct
                                if(i ~= length(app.data_url))
                                    while(s.bytes ~= 52428800)
                                        outfilename = websave(filename,url,options);
                                        s = dir(filename);
                                        pause(1);
                                    end
                                end
                                
                                % ------ Updating the downloaded file number in properties file
                                pause(1);
                                app.count = i;
                                change_xml_parameter(file_path,root_tab,[parameter_name],[string(i)],cell(0,0));
                            end
                            change_xml_parameter(file_path,root_tab,[parameter_name],["end"],cell(0,0));
                        catch
                            delete(f);
                            error_file =app.count + 1;
                            error_msg = srtcat("Download error in file: " , string(error_file));
                            errordlg(char(error_msg),'Error');
                            return;
                        end
                        pause(1);
                        jObj.setBusyText('Unpacking test data...');
                        try
                            % ------  Umpacking multi-pack Zip Package
                            source_7z = ['C:',filesep,'Program Files (x86)',filesep,'7-Zip',filesep,'7z.exe'];
                            surce_zip = [pwd,filesep,'H_HGGM_test_data.zip'];
                            
                            if(isfile(source_7z))
                                [status,result] = system(['"',source_7z,'"  x -y '  '"',surce_zip,'"' ]);
                                disp([status,result] );
                            else
                                % --- Downloading 7-Zip to install
                                zip_intaller_url =  'https://drive.google.com/open?id=1DJbV6_155hAvaTzmRyiewEPRXgUADWfI';
                                zip_installer_name = '7z1900.zip';
                                
                                outfilename = websave(zip_installer_name,zip_intaller_url);
                                pause(1);
                                exampleFiles = unzip(zip_installer_name,pwd);
                                pause(1);
                                % --- Installing 7-Zip in the PC
                                command = [pwd,filesep,'7z1900.exe']; %// external program; full path
                                [status,cmdout] = system([command ' & '])
                                
                                while(~isfile(source_7z))
                                    pause(2);
                                end
                                [status,result] = system(['"',source_7z,'"  x -y '  '"',surce_zip,'"' ]);
                                disp([status,result] );
                            end
                            
                        catch
                            delete(f);
                            errordlg('Unpackage error!!!','Error');
                            return;
                        end
                        jObj.stop;
                        jObj.setBusyText('All done!');
                        pause(2);
                        delete(f);
                        msgbox('Completed download!!!','Info');
                        result = true;
                    case 'Cancel'
                        result = false;
                        return;
                end
            else
                if(~simuling)
                    answer = questdlg('The test data has already been downloaded previously! Do you want to download again?', ...
                        'Select data', ...
                        'Yes I want','No','Cancel');
                    switch answer
                        case 'Yes I want'
                            change_xml_parameter(file_path,root_tab,[parameter_name],["1"],cell(0,0));
                            get_test_data(app,false);
                        case 'No'
                            result = true;
                            return;
                    end
                else
                    result = true;
                end
            end
        end
        
        function results = difine_paths(app)
            clc;
            addpath('common_functions');
            addpath('properties');
            addpath('EEG_ECOG');
            addpath('simulations/Sim1_hggm_convergence&jankova_conditions');
            addpath('simulations/Sim2_h_hggm_simplified_head_model');
            addpath('simulations/Sim3_h_hggm_realistic_head_model');
            addpath('simulations/Sim4_h_head_model_comparison');
            addpath('simulations/Sim_data');
            addpath('ssvep');
        end
    end
    

    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            clc;
            clear all;
            close all;
            addpath('common_functions');
            addpath('properties');
            addpath('EEG_ECOG');
            addpath('simulations/Sim1_hggm_convergence&jankova_conditions');
            addpath('simulations/Sim2_h_hggm_simplified_head_model');
            addpath('simulations/Sim3_h_hggm_realistic_head_model');
            addpath('simulations/Sim4_h_head_model_comparison');
            addpath('simulations/Sim_data');
            addpath('ssvep');
            
            
        end

        % Menu selected function: ExitMenu
        function ExitMenuSelected(app, event)
            delete(app) ;
        end

        % Menu selected function: HGGMConvergenceJankovaConditionsMenu
        function HGGMConvergenceJankovaConditionsMenuSelected(app, event)
            difine_paths(app)
            if(get_test_data(app,true))
                if(isempty( app.output_sourse) |  app.output_sourse==0)
                    app.output_sourse = uigetdir('tittle','Select the Output Folder');
                    if(app.output_sourse==0)
                        return;
                    end
                end
                Main_Jankova_test(app.output_sourse);
                msgbox('Completed operation','Info');
            end
        end

        % Menu selected function: SSVEPMenu
        function SSVEPMenuSelected(app, event)
            difine_paths(app);
            if(get_test_data(app,true))
                if(isempty( app.output_sourse) |  app.output_sourse==0)
                    app.output_sourse = uigetdir('tittle','Select the Output Folder');
                    if(app.output_sourse==0)
                        return;
                    end
                end
                Main_SSVEP_AnalyzeData(app.output_sourse);
                msgbox('Completed operation','Info');
            end
        end

        % Menu selected function: HHGGMSimplifiedHeadModelMenu
        function HHGGMSimplifiedHeadModelMenuSelected(app, event)
            difine_paths(app);
            if(get_test_data(app,true))
                if(isempty( app.output_sourse) |  app.output_sourse==0)
                    app.output_sourse = uigetdir('tittle','Select the Output Folder');
                    if(app.output_sourse==0)
                        return;
                    end
                end
                Main_simplified_em_penalty_test(app.output_sourse);
                msgbox('Completed operation','Info');
            end
        end

        % Menu selected function: HHGGMRealisticHeadModelMenu
        function HHGGMRealisticHeadModelMenuSelected(app, event)
            difine_paths(app);
            if(get_test_data(app,true))
                if(isempty( app.output_sourse) |  app.output_sourse==0)
                    app.output_sourse = uigetdir('tittle','Select the Output Folder');
                    if(app.output_sourse==0)
                        return;
                    end
                end
                Main_realistic_em_penalty_test(app.output_sourse);
                msgbox('Completed operation','Info');
            end
        end

        % Menu selected function: HHGGMHeadModelComparisonMenu
        function HHGGMHeadModelComparisonMenuSelected(app, event)
            difine_paths(app);
            if(get_test_data(app,true))
                if(isempty( app.output_sourse) |  app.output_sourse==0)
                    app.output_sourse = uigetdir('tittle','Select the Output Folder');
                    if(app.output_sourse==0)
                        return;
                    end
                end
                Main_H_HGGM_Head_Model_Comparison(app.output_sourse);
                msgbox('Completed operation','Info');
            end
        end

        % Menu selected function: SelectoutputfolderMenu
        function SelectoutputfolderMenuSelected(app, event)
            app.output_sourse = uigetdir('tittle','Select the Source Folder');
            if(app.output_sourse==0)
                return;
            end
            % create_data_structure(folder);
        end

        % Menu selected function: AllAnalyzeMenu
        function AllAnalyzeMenuSelected(app, event)
            difine_paths(app)
            if(get_test_data(app,true))
                if(isempty( app.output_sourse) |  app.output_sourse==0)
                    app.output_sourse = uigetdir('tittle','Select the Output Folder');
                    if(app.output_sourse==0)
                        return;
                    end
                end
                Main_Jankova_test(app.output_sourse);
                Main_simplified_em_penalty_test(app.output_sourse);
                Main_realistic_em_penalty_test(app.output_sourse);
                Main_H_HGGM_Head_Model_Comparison(app.output_sourse);
                Main_SSVEP_AnalyzeData(app.output_sourse);
                Main_EEG_ECOG_AnalyzeData(app.output_sourse);
                msgbox('Completed operation','Info');
            end
        end

        % Menu selected function: EEGECOGMenu
        function EEGECOGMenuSelected(app, event)
            difine_paths(app);
            if(get_test_data(app,true))
                if(isempty( app.output_sourse) |  app.output_sourse==0)
                    app.output_sourse = uigetdir('tittle','Select the Output Folder');
                    if(app.output_sourse==0)
                        return;
                    end
                end
                Main_EEG_ECOG_AnalyzeData(app.output_sourse);
                msgbox('Completed operation','Info');
            end
        end

        % Menu selected function: SelectdatafileMenu
        function SelectdatafileMenuSelected(app, event)
            
            get_test_data(app,false);
            
        end
    end

    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create HHGGMSimulationPackUIFigure
            app.HHGGMSimulationPackUIFigure = uifigure;
            app.HHGGMSimulationPackUIFigure.Position = [100 100 588 410];
            app.HHGGMSimulationPackUIFigure.Name = 'H HGGM Simulation Pack';

            % Create FileMenu
            app.FileMenu = uimenu(app.HHGGMSimulationPackUIFigure);
            app.FileMenu.Text = 'File';

            % Create SelectdatafileMenu
            app.SelectdatafileMenu = uimenu(app.FileMenu);
            app.SelectdatafileMenu.MenuSelectedFcn = createCallbackFcn(app, @SelectdatafileMenuSelected, true);
            app.SelectdatafileMenu.Text = 'Select data file';

            % Create SelectoutputfolderMenu
            app.SelectoutputfolderMenu = uimenu(app.FileMenu);
            app.SelectoutputfolderMenu.MenuSelectedFcn = createCallbackFcn(app, @SelectoutputfolderMenuSelected, true);
            app.SelectoutputfolderMenu.Text = 'Select output folder';

            % Create ExitMenu
            app.ExitMenu = uimenu(app.FileMenu);
            app.ExitMenu.MenuSelectedFcn = createCallbackFcn(app, @ExitMenuSelected, true);
            app.ExitMenu.Text = 'Exit';

            % Create ToolsMenu
            app.ToolsMenu = uimenu(app.HHGGMSimulationPackUIFigure);
            app.ToolsMenu.Text = 'Tools';

            % Create HGGMMenu
            app.HGGMMenu = uimenu(app.ToolsMenu);
            app.HGGMMenu.Text = 'HGGM';

            % Create HGGMConvergenceJankovaConditionsMenu
            app.HGGMConvergenceJankovaConditionsMenu = uimenu(app.HGGMMenu);
            app.HGGMConvergenceJankovaConditionsMenu.MenuSelectedFcn = createCallbackFcn(app, @HGGMConvergenceJankovaConditionsMenuSelected, true);
            app.HGGMConvergenceJankovaConditionsMenu.Text = 'HGGM Convergence & Jankova Conditions';

            % Create HHGGMSimplifiedHeadModelMenu
            app.HHGGMSimplifiedHeadModelMenu = uimenu(app.HGGMMenu);
            app.HHGGMSimplifiedHeadModelMenu.MenuSelectedFcn = createCallbackFcn(app, @HHGGMSimplifiedHeadModelMenuSelected, true);
            app.HHGGMSimplifiedHeadModelMenu.Text = 'H-HGGM Simplified Head Model';

            % Create HHGGMRealisticHeadModelMenu
            app.HHGGMRealisticHeadModelMenu = uimenu(app.HGGMMenu);
            app.HHGGMRealisticHeadModelMenu.MenuSelectedFcn = createCallbackFcn(app, @HHGGMRealisticHeadModelMenuSelected, true);
            app.HHGGMRealisticHeadModelMenu.Text = 'H-HGGM Realistic Head Model';

            % Create HHGGMHeadModelComparisonMenu
            app.HHGGMHeadModelComparisonMenu = uimenu(app.HGGMMenu);
            app.HHGGMHeadModelComparisonMenu.MenuSelectedFcn = createCallbackFcn(app, @HHGGMHeadModelComparisonMenuSelected, true);
            app.HHGGMHeadModelComparisonMenu.Text = 'H-HGGM Head Model Comparison';

            % Create SSVEPMenu
            app.SSVEPMenu = uimenu(app.ToolsMenu);
            app.SSVEPMenu.MenuSelectedFcn = createCallbackFcn(app, @SSVEPMenuSelected, true);
            app.SSVEPMenu.Text = 'SSVEP';

            % Create EEGECOGMenu
            app.EEGECOGMenu = uimenu(app.ToolsMenu);
            app.EEGECOGMenu.MenuSelectedFcn = createCallbackFcn(app, @EEGECOGMenuSelected, true);
            app.EEGECOGMenu.Text = 'EEG-ECOG';

            % Create AllAnalyzeMenu
            app.AllAnalyzeMenu = uimenu(app.ToolsMenu);
            app.AllAnalyzeMenu.MenuSelectedFcn = createCallbackFcn(app, @AllAnalyzeMenuSelected, true);
            app.AllAnalyzeMenu.Text = 'All Analyze';

            % Create HelpMenu
            app.HelpMenu = uimenu(app.HHGGMSimulationPackUIFigure);
            app.HelpMenu.Text = 'Help';
        end
    end

    methods (Access = public)

        % Construct app
        function app = h_hggm_simpack

            % Create and configure components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.HHGGMSimulationPackUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.HHGGMSimulationPackUIFigure)
        end
    end
end