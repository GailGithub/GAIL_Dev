function publishMathJax(filename,varargin)
% publish MATLAB scripts using the XSL file supporting MathJax
if nargin >= 3
   formatCompact = varargin{2};
else
   formatCompact = true;
   if nargin >= 2
      opts=varargin{1};
   end
end   
formatNow = get(0,'FormatSpacing'); %save existing format
if formatCompact
   format compact
else
   format loose
end
opts.stylesheet = 'mxdom2mathjaxbigfont.xsl';
publish(filename,opts);
set(0,'FormatSpacing',formatNow); %restore existing format

