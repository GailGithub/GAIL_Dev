function makeMovieRotateCoord(xPts,vidFileName,coord,nFrame,xDim,plotStyle)
% MAKEMOVIEROTATECOORD rotates the 3-D plot axes through several angles and
% outputs a movie

%% Parse inputs and set defaults
if nargin < 1
   error('Must have at least one input')
else
   d = size(xPts,2);
end
if nargin < 6
   plotStyle = '';
end
if nargin < 5
   xDim = [zeros(1,d); ones(1,d)];
elseif size(xDim,1) ~= 2
   xDim = [zeros(1,d); ones(1,d)];
end
if nargin < 4
   nFrame = 20; %number of frames per coordinate triple
elseif numel(nFrame) == 0
   nFrame = 20; %number of frames per coordinate triple
end
if nargin < 3
   coord = 1:d; %default way to go through the coordinates
elseif numel(coord) == 0
   coord = 1:d; %default way to go through the coordinates
end   
if nargin < 2
   vidFileName = 'rotatePts.mp4';
elseif ~ischar(vidFileName)
   vidFileName = 'rotatePts.mp4';
end

%% Initialize variables and display
gail.InitializeDisplay %initialze the display
frameInc = 0:1/(nFrame-1):1;       %frame incrementing
azInc = 90*frameInc;               %aximuth increments
elInc = 90*frameInc(end:-1:1);     %elevation increments
nLoop = length(coord) - 2;         %number of loops
figure;                            %initialzie figure
position = get(gcf,'Position');    %get position of figure
position(3) = position(4);         %and make it
set(gcf,'Position',position);      %square
set(gcf,'Color',[1 1 1]);          %make background color white
mov(nLoop*nFrame) = getframe(gcf); %initialize movie
k = 0;                             %initialize frame index

%% Creater movie frames
for l = 1:nLoop
   cd = coord(l:l+2);              %get coordianate triple
   h = plot3(xPts(:,cd(1)),xPts(:,cd(2)),xPts(:,cd(3)),'.'); %plot points
   if numel(plotStyle)
      eval(['set(h,' plotStyle ')']);
   end
   xlabel(['\(x_' int2str(cd(1)) '\)']) %label
   ylabel(['\(x_' int2str(cd(2)) '\)']) %all three
   zlabel(['\(x_' int2str(cd(3)) '\)']) %axes
   axis(xDim(:))                   %make axes have correct dimension
   axis('square')                  %make axes square
   set(gca,'Xtick',[xDim(1) (xDim(1)+xDim(2))/2 xDim(2)], ...
      'Ytick',[xDim(3) (xDim(3)+xDim(4))/2 xDim(4)], ...
      'Ztick',[xDim(5) (xDim(5)+xDim(6))/2 xDim(6)])
   for i = 1:nFrame
      k = k+1;                     %increment frame index
      view(azInc(i),elInc(i));     %set view
      mov(k) = getframe(gcf);      %record movie frame
   end
end

%% Write video file
vid = VideoWriter(vidFileName,'MPEG-4');
vid.FrameRate = 4;
open(vid)
writeVideo(vid,mov)
close(vid)

end