function res = Make_X_VA_Damp_Plot(trials, targetDirNum, dampNum, varargin)
    
    % Input Parser to add 'ShowViolations' option
    p = inputParser;
    argName = 'ShowViolations';
    defaultVal = false;
    addParameter(p, argName, defaultVal);

    parse(p,varargin{:});
    showViolations = p.Results.ShowViolations;

    data = FilterAndAnalyzeData(trials, targetDirNum, dampNum);
    
    % Line Color RGB vals
    lightGreyColor = [179, 179, 179]/256;
    lightGreenColor = [66, 245, 69]/256;
    
    % Plot
    fig1 = figure;
    set(gcf,'Color',[1,1,1]);
    
    % Position
    subplot(3,1,1);
    plot(data.x_data_wov, 'Color', lightGreyColor)
    hold on
    if (showViolations)
        plot(data.x_data_wv, 'Color', lightGreenColor)
    end
    plot(data.x_data_wov_mean, 'Color', 'k')
    plot(data.x_data_wov_ub, 'Color', 'k', 'LineStyle', '--');
    plot(data.x_data_wov_lb, 'Color', 'k', 'LineStyle', '--');
    hold off
    ylabel('x [m]', 'Interpreter', 'latex');
    title(sprintf('Damping: %s        Direction: %s', ...
                    data.DampingText, ...
                    data.TargetDirText));
	xlim([0, data.nSamples]);
    box('off');
                
    % VA
    subplot(3,1,2);
    plot(data.va_data_wov, 'Color', lightGreyColor)
    hold on 
    if (showViolations)
        plot(data.va_data_wv, 'Color', lightGreenColor)
    end
    plot(data.va_data_wov_mean, 'Color', 'k')
    plot(data.va_data_wov_ub, 'Color', 'k', 'LineStyle', '--');
    plot(data.va_data_wov_lb, 'Color', 'k', 'LineStyle', '--');
    hold off
    ylabel('$\dot{x} \ddot{x}$', 'Interpreter', 'latex');
    xlim([0, data.nSamples]);
    box('off');

    % Damping
    subplot(3,1,3);
    plot(data.damp_data_wov, 'Color', lightGreyColor)
    hold on
    if (showViolations)
        plot(data.damp_data_wv, 'Color', lightGreenColor)
    end
    plot(data.damp_data_wov_mean, 'Color', 'k')
    plot(data.damp_data_wov_ub, 'Color', 'k', 'LineStyle', '--');
    plot(data.damp_data_wov_lb, 'Color', 'k', 'LineStyle', '--');
    hold off
    ylabel('Damping [Ns/m]', 'Interpreter', 'latex');
    xlim([0, data.nSamples]);
    box('off');
    
    res = [fig1];
end