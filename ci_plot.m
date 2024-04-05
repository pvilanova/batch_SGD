function ciplot(lower,upper,x,colour)
% Plots a shaded region on a graph between specified lower and upper confidence intervals (L and U).

if length(lower)~=length(upper)
    error('lower and upper vectors must be same length')
end

% convert to row vectors so fliplr can work
if find(size(x)==(max(size(x))))<2
x=x'; end
if find(size(lower)==(max(size(lower))))<2
lower=lower'; end
if find(size(upper)==(max(size(upper))))<2
upper=upper'; end


h = fill([x fliplr(x)],[upper fliplr(lower)],colour,'EdgeColor','none','FaceAlpha',0.5,'DisplayName','CI');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

% hleg = legend('show');  % remove legend entry
% hleg.String(end) = [];
% legend('off');
