function figure_pdf_cropped(current_fig,nom,only)
% 
% if nargin==3, warning('third argument not used'), end
% disp('printing pdf...')
% save2pdf(nom,current_fig,600)
% disp('ok')
% 
% return


%set(gca,'FontSize',20,'FontName','Times');
%axis tight

if nargin==2
    set(current_fig,'PaperPositionMode','auto')
    disp('save fig...')
    saveas(current_fig,nom,'fig')
    
    disp('save tif...')
    saveas(current_fig,nom,'tif')
    
    disp('save pdf...') 
    set(current_fig,'Units','Inches');
    pos = get(current_fig,'Position');
    set(current_fig,'PaperPositionMode','Auto', ...
        'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(current_fig,nom,'-dpdf','-r300')
       
%     disp('save eps...')    
%     print(current_fig,'-depsc2',nom)
%     
%     disp('save pdf...')  
%     assert(not(system(sprintf('epstopdf %s.eps',nom))))
%     
%     disp(sprintf('delete %s.eps',nom))
%     delete ([nom '.eps'])
    
    %print(current_fig,'-dpdf',nom)
    
else
    set(current_fig,'PaperPositionMode','auto')
    if strcmp(only,'fig')
        disp('save fig...')
        saveas(current_fig,nom,'fig')
    elseif strcmp(only,'png')
        disp('save png...')
        saveas(current_fig,nom,'png')
    elseif strcmp(only,'pdf')
        disp('save pdf...')
         print(current_fig,nom,'-dpdf','-r300')
        %print(current_fig,'-depsc2','-r75',nom)
        %print(current_fig,'-depsc','-tiff','-r300',nom)
        %print('-depsc -tiff -r150', nom)
        %assert(not(system(sprintf('epstopdf %s.eps',nom))))
        %delete([nom '.eps']);
    end
end


