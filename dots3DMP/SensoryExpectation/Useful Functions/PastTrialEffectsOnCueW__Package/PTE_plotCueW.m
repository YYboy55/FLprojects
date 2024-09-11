%% Current functional function
function ax = PTE_plotCueW(vesW, vesW_avg, markerShape, clr,edgClr, xTickLabels)
    numColumns = size(vesW_avg{1}, 2); % num cohs
    numRows = size(vesW{1}, 1); % num bootstraps
    numGroups = size(xTickLabels,2); % Diff groups will be plotted in spaced out clusters along x-axis. MUST equal size(vesW,1) aka. num rows in the vesW cell array
    
    figure();
    hold on;
    for g = 1:numGroups
        for c = 1:length(vesW_avg) % Count num of cells i.e. num of trial types
            inGroupSpacing = c/800;
            % for plotting same data twice in same location on plot, (used to create star symbol from 2 triangles..)
%             if c==14
%                 inGroupSpacing = 7/295;
%             elseif c==13
%                 inGroupSpacing = 6/295;
%             elseif c==12
%                 inGroupSpacing = 5/295;
%             end
            if g==1
            % Plot data from vesW
            for col = 1:numColumns % num cohs of in each trial type (num columns in each cell)
%                 x = g/40 * col * ones(numRows, 1) + inGroupSpacing;
%                 y = vesW{g,c}(:,col);
%                 plot(x, y, '.', 'MarkerSize', 6, 'MarkerEdgeColor', clr{g,c});
            % Plot data from vesW_avg
                x = g/1 * col * ones(1, 1) + inGroupSpacing;
                y = vesW_avg{g,c}(:,col);
                plot(x, y, 'Marker',markerShape{g,c}, 'MarkerSize', 11, 'MarkerEdgeColor', edgClr{g,c},'MarkerFaceColor',clr{g,c});
            end
            else % Not group 1
                for col = 1:numColumns % num cohs of in each trial type (num columns in each cell)
%                     x = 1/40*(col)*ones(numRows, 1) + x1(g-1)/5 + inGroupSpacing; % Space between groups is x2 that of within group spacing
%                     y = vesW{g,c}(:,col);
%                     plot(x, y, '.', 'MarkerSize', 6, 'MarkerEdgeColor', clr{g,c});
                    % Plot data from vesW_avg
                    x = 1/1*(col)*ones(1, 1) + x1(g-1)/5 + inGroupSpacing;
                    y = vesW_avg{g,c}(:,col);
                    plot(x, y, 'Marker',markerShape{g,c}, 'MarkerSize', 11, 'MarkerEdgeColor', edgClr{g,c},'MarkerFaceColor',clr{g,c});
                end
            end
            if c == 1
                x1(g) = x;
            elseif c==length(vesW_avg)
                xLast = x;
            end
%             if g==1 && c==length(vesW)
%                 inGroupSpacing = (xLast-x1);
%             end
        end
        Ticks(g) = (xLast+x1(g))/2; % g/20 + ((c-1)/2)/200; %  ax.XLim = [mean(ax.XTick,2)/2,mean(ax.XTick,2)+mean(ax.XTick,2)/2];   % [.75,1.5];
    end
    ax = gca;
    ax.XTick = Ticks;
    ax.XLim = [Ticks(1)-inGroupSpacing Ticks(g)+inGroupSpacing]; % [mean(Ticks)/2,mean(Ticks)+mean(Ticks)/2];
    ax.YLim = [-.2,1.2];
    ax.XTickLabels = xTickLabels;

    title('Trial history affects cue weighting behavior'); ylabel('Vestibular cue weight');
end


%% Below should also be functional and almost identical, leaving here in case above doesn't work well :~]
% function ax = PTE_plotCueW(vesW, vesW_avg, clr, xTickLabels)
%     numColumns = size(vesW{1}, 2); % num cohs --> Attempting redesign where length(xTickLabels) matches num of diff cohs
%     numRows = size(vesW{1}, 1); % num bootstraps
%     numGroups = size(xTickLabels,2); % Should equal size(vesW,1) aka. num rows in the vesW cell array
%     figure();
%     hold on;
%     for g = 1:numGroups
%         for c = 1:length(vesW) % Count num of cells
%             % Plot data from vesW
%             for col = 1:numColumns
%                 x = g/5 * col * ones(numRows, 1) + c/100 - 1/100;
%                 y = vesW{g,c}(:,col);
%                 plot(x, y, '.', 'MarkerSize', 6, 'MarkerEdgeColor', clr{g,c});
%             % Plot data from vesW_avg
%                 x = g/5 * col * ones(1, 1) + c/100 - 1/100;
%                 y = vesW_avg{g,c}(:,col);
%                 plot(x, y, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', clr{g,c});
%             end
%         end
%         Ticks(g) = g/5 + ((c-1)/2)/100; %  ax.XLim = [mean(ax.XTick,2)/2,mean(ax.XTick,2)+mean(ax.XTick,2)/2];   % [.75,1.5];
%     end
%     ax = gca;
%     ax.XTick = Ticks;
%     ax.XLim = [mean(Ticks)/2,mean(Ticks)+mean(Ticks)/2];
%     ax.YLim = [-.1,1.1];
%     ax.XTickLabels = xTickLabels;
% 
%     title('Trial history affects cue weighting behavior'); ylabel('Vestibular cue weight');
% end

%% For plotting multiple cohs
% function ax = PTE_plotCueW(vesW, vesW_avg, clr, xTickLabels)
%     numColumns = size(vesW{1}, 2); % num cohs --> Attempting redesign where length(xTickLabels) matches num of diff cohs
%     numRows = size(vesW{1}, 1); % num bootstraps
%     for c = 1:length(vesW) % Count num of cells
%         % Plot data from vesW
%         for col = 1:numColumns
%             x = col * ones(numRows, 1) + c/50 - 1/50;
%             y = vesW{c}(:,col);
%             plot(x, y, '.', 'MarkerSize', 6, 'MarkerEdgeColor', clr{c});
%         end
% 
%         % Plot data from vesW_avg
%         for col = 1:numColumns
%             x = col * ones(1, 1) + c/50 - 1/50;
%             y = vesW_avg{c}(:,col);
%             plot(x, y, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', clr{c});
%         end
%     end
%     hold off;
%     ax = gca;
% 
%     for g = 1:xTickLabels
%         ax.XTick = g + ((c-1)/2)/50; %  ax.XLim = [mean(ax.XTick,2)/2,mean(ax.XTick,2)+mean(ax.XTick,2)/2];   % [.75,1.5];
%         ax.XLim = [.75,1.5];
%         ax.YLim = [-.1,1.1];
%         ax.XTickLabels = xTickLabels{g};
%     end
%     title('Trial history affects cue weighting behavior'); ylabel('Vestibular cue weight');
% end

%%
% 
% if numColumns == 1% Plotting for only a single coherence
%     ax.XTick = 1 + ((c-1)/2)/10;
%     ax.XLim = [.75,1.5];
%     ax.YLim = [-.25,1.3];
%     ax.XTickLabels = xTickLabels;
% elseif numColumns == 2 % For 2 coherences
%     ax.XTick = [1 + ((c-1)/2)/10, 2 + ((c-1)/2)/10];
%     ax.XTickLabels = xTickLabels;
%     ax.XLim = [.75,2.45];
%     ax.YLim = [-.25,1.3];
%     title('Trial history affects cue weighting behavior'); ylabel('Vestibular cue weight');
% end
