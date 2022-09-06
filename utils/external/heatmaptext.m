function  heatmaptext(data,varargin)
%HEATMAPTEXT creates a heatmap with the values shown as text.
%
%   HEATMAPTEXT(SEQ) counts the number of occurrences of each codon in the
%   sequence and displays a formatted table of the result.
%
%   HEATMAPTEXT(...,'PRECISION',P) displays the data values with a maximum N
%   digits of precision.  The default number of digits of precision is 2.
%
%   HEATMAPTEXT(SEQ, ... ,'FONTCOLOR',COL) sets the color of the text to
%   COL, where COL is a valied Handle Graphics color property.
%
%   HEATMAPTEXT(SEQ, ... ,'COLORBAR',false) hides the colorbar.
%
%   Examples:
%
%       ih = gallery('invhess',10);
%       heatmaptext(ih);
%
%       figure;
%       m = gallery('moler',30);
%       heatmaptext(m);
%
%       figure;
%       heatmaptext(rand(20,10),'fontcolor','r','precision',4);
%       colormap(bone);
%
%   See also COLORMAP, IMAGESC, TEXT.

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 2007/08/08 21:56:22 $

precision = 2;
textColor = 'k';
cBar = true;
if nargin > 1
    if rem(nargin,2)== 0
        error('heatmaptext:IncorrectNumberOfArguments',...
            'Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'precision','fontcolor','colorbar'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = strmatch(lower(pname), okargs);
        if isempty(k)
            error('heatmaptext:UnknownParameterName',...
                'Unknown parameter name: %s.',pname);
        elseif length(k)>1
            error('heatmaptext:AmbiguousParameterName',...
                'Ambiguous parameter name: %s.',pname);
        else
            switch(k)
                case 1  % frame
                    precision = pval;
                case 2  % direction
                    textColor = pval;
                case 3
                    cBar = pval;

            end
        end
    end
end


im = imagesc(data);
imAxes =get(im,'parent');
hFig = get(imAxes,'parent');
fs = getBestFontSize(imAxes);
    showText = true;
if fs == 0
    showText = false;
    fs = 9;
end
axis off;
if cBar
    colorbar;
end
[rows, cols] = size(data);
if max(rows,cols) > 50
    warning('This might be quite slow...');
end

textHandles = zeros(size(data))';
for i = 1:rows
    for j = 1:cols
                 textHandles(j,i) = text(j,i,num2str(data(i,j),precision),...
                'color',textColor,'horizontalAlignment','center',...
                'fontsize',fs,'clipping','on','visible','off');
    end
end
if  showText 
    set(textHandles,'visible','on');
end
set(imAxes,'UserData',textHandles);
localSetupPositionListener(hFig,imAxes);
localSetupLimitListener(imAxes);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fs = getBestFontSize(imAxes)
% Try to keep font size reasonable for text
hFig = get(imAxes,'Parent');
magicNumber = 80;
nrows = diff(get(imAxes,'YLim'));
ncols = diff(get(imAxes,'XLim'));
if ncols < magicNumber && nrows < magicNumber
    ratio = max(get(hFig,'Position').*[0 0 0 1])/max(nrows,ncols);
elseif ncols < magicNumber
    ratio = max(get(hFig,'Position').*[0 0 0 1])/ncols;
elseif nrows < magicNumber
    ratio = max(get(hFig,'Position').*[0 0 0 1])/nrows;
else
    ratio = 1;
end
fs = min(9,ceil(ratio/4));    % the gold formula
if fs < 4
    fs = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localSetupPositionListener(hFig,imAxes)
% helper function to sets up listeners for resizing, so we can detect if
% we would need to change the fontsize
hgp     = findpackage('hg');
figC   = findclass(hgp,'figure');
PostPositionListener = handle.listener(hFig,'ResizeEvent',...
    {@localPostPositionListener,imAxes});
% PostPositionListener = handle.listener(hFig,figC.findprop('Position'),...
%     'PropertyPostSet',{@localPostPositionListener,imAxes});
% 

listeners = getappdata(hFig,'HeatMapListeners');

setappdata(hFig,'HeatMapListeners',[listeners PostPositionListener]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localPostPositionListener(hSrc,event,imAxes) %#ok
textHandles = get(imAxes,'UserData');
fs = getBestFontSize(imAxes);
if fs > 0
    set(textHandles,'fontsize',fs,'visible','on');
else
    set(textHandles,'visible','off');
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localSetupLimitListener(imAxes)
% helper function to sets up listeners for zooming, so we can detect if
% we would need to change the fontsize
hgp     = findpackage('hg');
axesC   = findclass(hgp,'axes');
LimListener = handle.listener(imAxes,[axesC.findprop('XLim') axesC.findprop('YLim')],...
    'PropertyPostSet',@localLimitListener);

hFig = get(imAxes,'Parent');
listeners = getappdata(hFig,'HeatMapListeners');

setappdata(hFig,'HeatMapListeners',[listeners LimListener]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localLimitListener(hSrc,event)%#ok

imAxes = event.AffectedObject;

textHandles = get(imAxes,'UserData');
fs = getBestFontSize(imAxes);
if fs > 0
    set(textHandles,'fontsize',fs,'visible','on');
else
    set(textHandles,'visible','off');
end

