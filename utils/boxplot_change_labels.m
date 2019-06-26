function boxplot_change_labels(final_res, tickLabelStr, fontSize)

nn=length(tickLabelStr);

    % group boxes
    width = 0.5;
    sh = 0.1; 
    
    pos = [];
    for i=1:nn
        pos = [pos i+sh];
    end
    
    wid = width * ones(1,length(pos));
    
    % boxplot
    %figure
    boxplot(final_res, ...%'notch', 'on', ...
            'positions', pos,...
            'widths', wid) 
    
    % label, change fontsize
    text_h = findobj(gca, 'Type', 'text');
    rotation = 45;
    
    for cnt = 1:length(text_h)
        set(text_h(cnt),    'FontSize', fontSize,...
                            'Rotation', rotation, ...
                            'String', tickLabelStr{length(tickLabelStr)-cnt+1}, ...
                            'HorizontalAlignment', 'right')
    end
    
    % 'VerticalAlignment', 'cap', ...
    
    % smaller box for axes, in order to un-hide the labels
    squeeze = 0.2;
    left = 0.02;
    right = 1;
    bottom = squeeze;
    top = 1-squeeze;
    
    %set(gca, 'OuterPosition', [left bottom right top])
        
    % remove outliers
    hout = findobj(gca,'tag','Outliers');
    for out_cnt = 1 : length(hout)
        set(hout(out_cnt), 'Visible', 'off')
    end