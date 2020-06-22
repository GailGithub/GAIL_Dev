function colorMarkerLinePlot(h,dd,colorSequence,markerSequence,markerSize, ...
   lineSequence)
nColor = length(colorSequence);
nMarker = length(markerSequence);
nLine = length(lineSequence);
whColor = colorSequence{mod(dd-1,nColor)+1};
whMarker = markerSequence{mod(dd-1,nMarker)+1};
whSize = markerSize{mod(dd-1,nMarker)+1};
whLine = lineSequence{mod(dd-1,nLine)+1};
set(h,'color',whColor,'marker',whMarker,'MarkerFaceColor',whColor, ...
   'MarkerSize',whSize,'LineStyle',whLine)
