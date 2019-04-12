%%
% Authors:
% - Deirel Paz Linares
% - Eduardo Gonzalez Moreira
% - Pedro A. Valdes Sosa

% Date: March 16, 2019



% Updates
% - Ariosky Areces Gonzalez


classdef h_hggm_simpack < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        HHGGMSimulationPackUIFigure   matlab.ui.Figure
        FileMenu                      matlab.ui.container.Menu
        OutputsourseMenu              matlab.ui.container.Menu
        ExitMenu                      matlab.ui.container.Menu
        ToolsMenu                     matlab.ui.container.Menu
        HGGMMenu                      matlab.ui.container.Menu
        HGGMConvergenceJankovaConditionsMenu  matlab.ui.container.Menu
        HHGGMSimplifiedHeadModelMenu  matlab.ui.container.Menu
        HHGGMRealisticHeadModelMenu   matlab.ui.container.Menu
        HHGGMHeadModelComparisonMenu  matlab.ui.container.Menu
        SSVEPMenu                     matlab.ui.container.Menu
        AllAnalyzeMenu                matlab.ui.container.Menu
        HelpMenu                      matlab.ui.container.Menu
    end

    
    properties (Access = public)
        output_sourse % Description
    end
    

    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            clc;
            clear all;
            close all;
            addpath('simulations/Sim_Common_Functions');
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
            if(isempty( app.output_sourse) |  app.output_sourse==0)
                app.output_sourse = uigetdir('tittle','Select the Output Folder');
                if(app.output_sourse==0)
                    return;
                end
            end
            Main_Jankova_test(app.output_sourse);
            msgbox('Completed operation','Info');
        end

        % Menu selected function: SSVEPMenu
        function SSVEPMenuSelected(app, event)
           if(isempty( app.output_sourse) |  app.output_sourse==0)
                app.output_sourse = uigetdir('tittle','Select the Output Folder');
                if(app.output_sourse==0)
                    return;
                end
            end
            Main_SSVEP_AnalyzeData(app.output_sourse);
            msgbox('Completed operation','Info');
        end

        % Menu selected function: HHGGMSimplifiedHeadModelMenu
        function HHGGMSimplifiedHeadModelMenuSelected(app, event)
           if(isempty( app.output_sourse) |  app.output_sourse==0)
                app.output_sourse = uigetdir('tittle','Select the Output Folder');
                if(app.output_sourse==0)
                    return;
                end
            end
            Main_simplified_em_penalty_test(app.output_sourse);
            msgbox('Completed operation','Info');
        end

        % Menu selected function: HHGGMRealisticHeadModelMenu
        function HHGGMRealisticHeadModelMenuSelected(app, event)
           if(isempty( app.output_sourse) |  app.output_sourse==0)
                app.output_sourse = uigetdir('tittle','Select the Output Folder');
                if(app.output_sourse==0)
                    return;
                end
            end
            Main_realistic_em_penalty_test(app.output_sourse);
            msgbox('Completed operation','Info');
        end

        % Menu selected function: HHGGMHeadModelComparisonMenu
        function HHGGMHeadModelComparisonMenuSelected(app, event)
           if(isempty( app.output_sourse) |  app.output_sourse==0)
                app.output_sourse = uigetdir('tittle','Select the Output Folder');
                if(app.output_sourse==0)
                    return;
                end
            end
            Main_H_HGGM_Head_Model_Comparison(app.output_sourse);
            msgbox('Completed operation','Info');
        end

        % Menu selected function: OutputsourseMenu
        function OutputsourseMenuSelected(app, event)
            app.output_sourse = uigetdir('tittle','Select the Source Folder');
            if(app.output_sourse==0)
                return;
            end
            % create_data_structure(folder);
        end

        % Menu selected function: AllAnalyzeMenu
        function AllAnalyzeMenuSelected(app, event)
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
            msgbox('Completed operation','Info');
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

            % Create OutputsourseMenu
            app.OutputsourseMenu = uimenu(app.FileMenu);
            app.OutputsourseMenu.MenuSelectedFcn = createCallbackFcn(app, @OutputsourseMenuSelected, true);
            app.OutputsourseMenu.Text = 'Output sourse';

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