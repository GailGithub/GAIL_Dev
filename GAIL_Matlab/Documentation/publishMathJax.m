function publishMathJax(filename,varargin)
% publish MATLAB scripts using the XSL file supporting MathJax
if nargin >= 3
   formatcompact = varargin{2};
else
   formatcompact = true;
   if nargin >= 2
      opts=varargin{1};
   end
end   
if formatcompact
   format compact
end
opts.stylesheet = 'mxdom2mathjaxbigfont.xsl';
publish(filename,opts);

