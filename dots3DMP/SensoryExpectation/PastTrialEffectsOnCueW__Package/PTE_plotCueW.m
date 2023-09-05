function ax = PTE_plotCueW(vesW, vesW_avg, clr)
    numColumns = size(vesW{1}, 2); % num cohs
    numRows = size(vesW{1}, 1); % num bootstraps

    figure;
    hold on;
    for c = 1:length(vesW) % Count num of cells
        % Plot data from vesW
        for col = 1:numColumns
            x = col * ones(numRows, 1) + c/10 - 1/10;
            y = vesW{c}(:,col);
            plot(x, y, '.', 'MarkerSize', 6, 'MarkerEdgeColor', clr{c});
        end

        % Plot data from vesW_avg
        for col = 1:numColumns
            x = col * ones(1, 1) + c/10 - 1/10;
            y = vesW_avg{c}(:,col);
            plot(x, y, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', clr{c});
        end
    end
    hold off;
    ax = gca;
    
    % Plotting for only a single coherence
%     ax.XTick = 1 + ((c-1)/2)/10;
%     ax.XLim = [.75,1.5];
%     ax.XTickLabels = {'High coh trials'};
% For 2 coherences
    ax.XTick = [1 + ((c-1)/2)/10, 2 + ((c-1)/2)/10];
    ax.XTickLabels = {'Low coh trials','High coh trials'};
    ax.XLim = [.75,2.45];
    ax.YLim = [-.25,1.3];
    title('Trial history affects cue weighting behavior'); ylabel('Vestibular cue weight');
end