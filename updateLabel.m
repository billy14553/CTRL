function  updateLabel(alllabel,dir)
% update label
%dir 0----xlabel,1----ylabel
if dir==0
    mytick = @xticks;
    myticklabel = @xticklabels;
else
    mytick = @yticks;
    myticklabel = @yticklabels;
end
%updateLabel
if ~isempty(alllabel)
    xtk = mytick();
    xlabels =myticklabel();
    for i = 1:1:length(xtk)
        if or(xtk(i)<1,xtk(i)>length(alllabel))
              xlabels{i} = "";
        else
             xlabels{i} = num2str(alllabel(xtk(i)));
        end
    end
    myticklabel(xlabels);
end
end

