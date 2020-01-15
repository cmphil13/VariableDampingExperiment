function res = Make_X_VA_Damp_Plot_Combined(trials, targetDirNum, varargin)
    
    % Input Parser to add 'ShowViolations' option
    p = inputParser;
    argName = 'ShowViolations';
    defaultVal = false;
    addParameter(p, argName, defaultVal);
    addParameter(p, 'ShowMovementMatch', false);
    
    
    parse(p,varargin{:});
    showViolations = p.Results.ShowViolations;
    showMovementMatch = p.Results.ShowMovementMatch;
    
    targetDirText = {'Left', 'Right', 'Backward', 'Forward'};

    % Plot
    fig1 = figure;
    set(fig1,'Color',[1,1,1]);
    
    if (showMovementMatch)
        fig2 = figure;
        set(fig2, 'Color', [1,1,1]);
    end
    
    orange = [255, 165, 0]/256;
    dampColors = {'b', 'g', orange};
    dampTexts = {};
    
    % Make Plots
    for dampNum = 1:3

        data = FilterAndAnalyzeData(trials, targetDirNum, dampNum);
    
        % Line Color RGB vals
        lightGreyColor = [179, 179, 179]/256;
        lightGreenColor = [66, 245, 69]/256;
        
        % set fig1 as current figure
        set(0, 'CurrentFigure', fig1);

        % Position
        subplot(3,4,1+(dampNum-1));
%         plot(data.wov.x.data_, 'Color', lightGreyColor)
%         hold on
        ub = data.wov.x.mean_ + data.wov.x.std_;
        lb = data.wov.x.mean_ - data.wov.x.std_;
        plot(data.wov.x.mean_, 'Color', 'k')
        hold on
        if (showViolations)
            plot(data.wv.x.data_, 'Color', lightGreenColor)
        end
        plot(ub, 'Color', 'k', 'LineStyle', '--');
        plot(lb, 'Color', 'k', 'LineStyle', '--');
        hold off
        if (dampNum == 1)
            ylabel('x [m]', 'Interpreter', 'latex');
        end
        if (dampNum == 2)
            title(sprintf('Direction: %s\n\nDamping: %s', ...
                        targetDirText{targetDirNum}, ...
                        data.DampingText));
        else
            title(sprintf('Damping: %s', ...
                          data.DampingText));
        end

        xlim([0, data.nSamples]);
        box('off');
                
        % VA
        subplot(3,4,5+(dampNum-1));
%         plot(data.wov.va.data_, 'Color', lightGreyColor)
        ub = data.wov.va.mean_ + data.wov.va.std_; % makes 3 std's
        lb = data.wov.va.mean_ - data.wov.va.std_; % makes 3 std's
        
        plot(data.wov.va.mean_, 'Color', 'k')
        hold on
        if (showViolations)
            plot(data.wv.va.data_, 'Color', lightGreenColor)
        end

        plot(ub, 'Color', 'k', 'LineStyle', '--');
        plot(lb, 'Color', 'k', 'LineStyle', '--');
        hold off
        if (dampNum == 1)
            ylabel('$\dot{x} \ddot{x}$', 'Interpreter', 'latex');
        end
        xlim([0, data.nSamples]);
        box('off');

        % Damping
        subplot(3,4,9+(dampNum-1));
%         plot(data.wov.damp.data_, 'Color', lightGreyColor)
        ub = data.wov.damp.mean_ + data.wov.damp.std_;
        lb = data.wov.damp.mean_ - data.wov.damp.std_;

        plot(data.wov.damp.mean_, 'Color', 'k')
        hold on
        if (showViolations)
            plot(data.wv.damp.data_, 'Color', lightGreenColor)
        end
        plot(ub, 'Color', 'k', 'LineStyle', '--');
        plot(lb, 'Color', 'k', 'LineStyle', '--');
        hold off
        if (dampNum == 1)
            ylabel('Damping [Ns/m]', 'Interpreter', 'latex');
        end
        xlabel('Time [ms]', 'Interpreter', 'latex');
        xlim([0, data.nSamples]);
        box('off');
        
        % Means - x
        subplot(3,4,4)
        if (dampNum > 1)
            hold on
        end
        plot(data.wov.x.mean_, 'Color', dampColors{dampNum},'LineWidth', 1.5);
        hold off
        box('off');

        % Means - va
        subplot(3,4,8)
        if (dampNum > 1)
            hold on
        end
        plot(data.wov.va.mean_, 'Color', dampColors{dampNum},'LineWidth', 1.5);
        hold off
        box('off');

        % Means - damp
        subplot(3,4,12)
        if (dampNum > 1)
            hold on
        end
        plot(data.wov.damp.mean_, 'Color', dampColors{dampNum},'LineWidth', 1.5);
        hold off
        box('off');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%            MOVEMENT MATCHED PLOT        %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (showMovementMatch)
            % set fig2 as current figure
            set(0, 'CurrentFigure', fig2);
            
            % Position
            subplot(3,4,1+(dampNum-1));
            plot(data.wov.x_mm.data_, 'Color', lightGreyColor)
            hold on
            if (showViolations)
                plot(data.wv.x_mm.data_, 'Color', lightGreenColor)
            end
            plot(data.wov.x_mm.mean_, 'Color', 'k')
            plot(data.wov.x_mm.ub_, 'Color', 'k', 'LineStyle', '--');
            plot(data.wov.x_mm.lb_, 'Color', 'k', 'LineStyle', '--');
            hold off
            if (dampNum == 1)
                ylabel('x [m]', 'Interpreter', 'latex');
            end
            if (dampNum == 2)
                title(sprintf('Direction: %s\n\nDamping: %s', ...
                            data.TargetDirText, ...
                            data.DampingText));
            else
                title(sprintf('Damping: %s', ...
                              data.DampingText));
            end

            xlim([0, data.nSamples]);
            box('off');

            % VA
            subplot(3,4,5+(dampNum-1));
            plot(data.wov.va_mm.data_, 'Color', lightGreyColor)
            hold on
            if (showViolations)
                plot(data.wv.va_mm.data_, 'Color', lightGreenColor)
            end
            plot(data.wov.va_mm.mean_, 'Color', 'k')
            plot(data.wov.va_mm.ub_, 'Color', 'k', 'LineStyle', '--');
            plot(data.wov.va_mm.lb_, 'Color', 'k', 'LineStyle', '--');
            hold off
            if (dampNum == 1)
                ylabel('$\dot{x} \ddot{x}$', 'Interpreter', 'latex');
            end
            xlim([0, data.nSamples]);
            box('off');

            % Damping
            subplot(3,4,9+(dampNum-1));
            plot(data.wov.damp_mm.data_, 'Color', lightGreyColor)
            hold on
            if (showViolations)
                plot(data.wv.damp_mm.data_, 'Color', lightGreenColor)
            end
            plot(data.wov.damp_mm.mean_, 'Color', 'k')
            plot(data.wov.damp_mm.ub_, 'Color', 'k', 'LineStyle', '--');
            plot(data.wov.damp_mm.lb_, 'Color', 'k', 'LineStyle', '--');
            hold off
            if (dampNum == 1)
                ylabel('Damping [Ns/m]', 'Interpreter', 'latex');
            end
            xlabel('Time [ms]', 'Interpreter', 'latex');
            xlim([0, data.nSamples]);
            box('off');

            % Means - x
            subplot(3,4,4)
            if (dampNum > 1)
                hold on
            end
            plot(data.wov.x_mm.mean_, dampColors{dampNum}, 'LineWidth', 1.5);
            hold off
            box('off');

            % Means - va
            subplot(3,4,8)
            if (dampNum > 1)
                hold on
            end
            plot(data.wov.va_mm.mean_, dampColors{dampNum}, 'LineWidth', 1.5);
            hold off
            box('off');

            % Means - damp
            subplot(3,4,12)
            if (dampNum > 1)
                hold on
            end
            plot(data.wov.damp_mm.mean_, dampColors{dampNum}, 'LineWidth', 1.5);
            hold off
            box('off');

            end

            % Add to dampText
            dampText{dampNum} = data.DampingText;
        
    end
    
    % Set fig1 as current figure
    set(0, 'CurrentFigure', fig1)
    
    % Apply legends
    subplot(3,4,4)
    legend(dampText,'FontSize',11)
    legend('boxoff');
    if (targetDirNum == 1 || targetDirNum == 3)
        legend('Location', 'southeast');
    else
        legend('Location', 'northeast');
    end
    
%     % Apply legends
%     subplot(3,4,8)
%     legend(dampText)
%     legend('boxoff');
%     if (targetDirNum == 1 || targetDirNum == 3)
%         legend('Location', 'southeast');
%     else
%         legend('Location', 'northeast');
%     end
% 
%     % Apply legends
%     subplot(3,4,12)
%     legend(dampText)
%     legend('boxoff');
%     if (targetDirNum == 1 || targetDirNum == 3)
%         legend('Location', 'southeast');
%     else
%         legend('Location', 'northeast');
%     end
    subplot(3,4,12)
    xlabel('Time [ms]', 'Interpreter', 'latex');
    
    
    % Set X - YLims Equal
    axes_arr = [subplot(3,4,1), ...
                subplot(3,4,2), ...
                subplot(3,4,3), ...
                subplot(3,4,4)];
    SetYLimsEqual(axes_arr, ...
                 'TopPadPercent', 0.1, ...
                 'BottomPadPercent', 0.1);

    % Set VA - YLims Equal
    axes_arr = [subplot(3,4,5), ...
                subplot(3,4,6), ...
                subplot(3,4,7), ...
                subplot(3,4,8)];
    SetYLimsEqual(axes_arr, ...
                 'TopPadPercent', 0.1, ...
                 'BottomPadPercent', 0.1);
    
    % Set Damping - YLims Equal
    axes_arr = [subplot(3,4,9), ...
                subplot(3,4,10), ...
                subplot(3,4,11), ...
                subplot(3,4,12)];
    SetYLimsEqual(axes_arr, ...
                 'TopPadPercent', 0.05, ...
                 'BottomPadPercent', 0.05);    
             
             
    if (showMovementMatch)
        % Set fig1 as current figure
        set(0, 'CurrentFigure', fig2)

        % Apply legends
        subplot(3,4,4)
        legend(dampText)
        legend('boxoff');
        if (targetDirNum == 1 || targetDirNum == 3)
            legend('Location', 'southeast');
        else
            legend('Location', 'northeast');
        end

        % Apply legends
        subplot(3,4,8)
        legend(dampText)
        legend('boxoff');
        if (targetDirNum == 1 || targetDirNum == 3)
            legend('Location', 'southeast');
        else
            legend('Location', 'northeast');
        end

        % Apply legends
        subplot(3,4,12)
        legend(dampText)
        legend('boxoff');
        if (targetDirNum == 1 || targetDirNum == 3)
            legend('Location', 'southeast');
        else
            legend('Location', 'northeast');
        end
        xlabel('Time [ms]', 'Interpreter', 'latex');


        % Set X - YLims Equal
        axes_arr = [subplot(3,4,1), ...
                    subplot(3,4,2), ...
                    subplot(3,4,3), ...
                    subplot(3,4,4)];
        SetYLimsEqual(axes_arr, ...
                     'TopPadPercent', 0.1, ...
                     'BottomPadPercent', 0.1);

        % Set VA - YLims Equal
        axes_arr = [subplot(3,4,5), ...
                    subplot(3,4,6), ...
                    subplot(3,4,7), ...
                    subplot(3,4,8)];
        SetYLimsEqual(axes_arr, ...
                     'TopPadPercent', 0.1, ...
                     'BottomPadPercent', 0.1);

        % Set Damping - YLims Equal
        axes_arr = [subplot(3,4,9), ...
                    subplot(3,4,10), ...
                    subplot(3,4,11), ...
                    subplot(3,4,12)];
        SetYLimsEqual(axes_arr, ...
                     'TopPadPercent', 0.05, ...
                     'BottomPadPercent', 0.05);    
    end
                 
                 
    res = [fig1];
end