classdef h_hggm_simpack < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        HHGGMSimulationPackUIFigure   matlab.ui.Figure
        FileMenu                      matlab.ui.container.Menu
        TestdatasourceMenu            matlab.ui.container.Menu
        OutputsourceMenu              matlab.ui.container.Menu
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
    end
    
    methods (Access = private)
        
        function result = get_test_data(app)
            result = false;
            if((isempty( app.test_data_sourse) |  app.test_data_sourse==0) & isempty(app.downloaded))
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
                                app.test_data_sourse = fullfile(path,file);
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
                            url = app.data_url;
                            filename = strcat('H_HGGM_test_data.zip');
                            options = weboptions('Timeout',Inf,'RequestMethod','get');
                            outfilename = websave(filename,url,options);
                        catch
                            delete(f);
                            errordlg('Download error!!!','Error');
                            return;
                        end
                        jObj.setBusyText('Unpacking test data...');
                        try
                            exampleFiles = unzip(filename,pwd);
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
                        app.downloaded = true;
                    case 'Cancel'
                        result = false;
                        return;
                end
            else
                result = true;
            end
        end
    end
    

    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            clc;
            clear all;
            close all;
            addpath('common_functions');
            addpath('EEG_ECOG');
            addpath('simulations/Sim1_hggm_convergence&jankova_conditions');
            addpath('simulations/Sim2_h_hggm_simplified_head_model');
            addpath('simulations/Sim3_h_hggm_realistic_head_model');
            addpath('simulations/Sim4_h_head_model_comparison');
            addpath('simulations/Sim_data');
            addpath('ssvep');
            
            app.downloaded = false;
            app.data_url = "https://drive.google.com/uc?id=1QLl5ZmCqVf9JDzhpxxUfFM6-4rmmtTM6";
        end

        % Menu selected function: ExitMenu
        function ExitMenuSelected(app, event)
            delete(app) ;
        end

        % Menu selected function: HGGMConvergenceJankovaConditionsMenu
        function HGGMConvergenceJankovaConditionsMenuSelected(app, event)
            if(get_test_data(app))
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
            if(get_test_data(app))
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
            if(get_test_data(app))
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
            if(get_test_data(app))
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
            if(get_test_data(app))
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

        % Menu selected function: OutputsourceMenu
        function OutputsourceMenuSelected(app, event)
            app.output_sourse = uigetdir('tittle','Select the Source Folder');
            if(app.output_sourse==0)
                return;
            end
            % create_data_structure(folder);
        end

        % Menu selected function: AllAnalyzeMenu
        function AllAnalyzeMenuSelected(app, event)
            if(get_test_data(app))
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
            if(get_test_data(app))
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

        % Menu selected function: TestdatasourceMenu
        function TestdatasourceMenuSelected(app, event)
            
            get_test_data(app);
            
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

            % Create TestdatasourceMenu
            app.TestdatasourceMenu = uimenu(app.FileMenu);
            app.TestdatasourceMenu.MenuSelectedFcn = createCallbackFcn(app, @TestdatasourceMenuSelected, true);
            app.TestdatasourceMenu.Text = 'Test data source';

            % Create OutputsourceMenu
            app.OutputsourceMenu = uimenu(app.FileMenu);
            app.OutputsourceMenu.MenuSelectedFcn = createCallbackFcn(app, @OutputsourceMenuSelected, true);
            app.OutputsourceMenu.Text = 'Output source';

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