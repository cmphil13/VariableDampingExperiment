function res = SetYLimsEqual(axes_array, varargin)
    
    % Input parser
    p = inputParser;
    addParameter(p, 'TopPadPercent', 0.05);
    addParameter(p, 'BottomPadPercent', 0.05);
    parse(p, varargin{:});
    
    tpp = p.Results.TopPadPercent;
    bpp = p.Results.BottomPadPercent;


    % Gather ylims
    ylim_arr = [];
    for i = 1:length(axes_array)
        ylim_arr = [ylim_arr; axes_array(i).YLim];
    end
    
    % Get min and max bounds
    ylim_bounds = [min(ylim_arr(:,1)), max(ylim_arr(:,2))];
    
    % Add extra space to bounds
    ylim_bounds_range = abs(ylim_bounds(2) - ylim_bounds(1));
    ylim_bounds(1) = ylim_bounds(1) - bpp*ylim_bounds_range;
    ylim_bounds(2) = ylim_bounds(2) + tpp*ylim_bounds_range;
    
    % Apply bounds to each axes
    for i = 1:length(axes_array)
        axes_array(i).YLim = ylim_bounds;
    end

end
