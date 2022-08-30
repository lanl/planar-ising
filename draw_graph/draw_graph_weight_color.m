% Jason K. Johnson
% LANL, Aug 2010
%
% SYNOPSIS draw an undirected graph with edge weights controlling
% shading of drawn edges.
% INPUTS:
% vx,vy: horizontal and vertical locations of vertices in the plane
% edges: m x 2 matrix of edges {edges(k,1),edges(k,2)}, k = 1,...,m
% w: edge weights
% cmap_nodes: node colors, n x 3 matrix of color codes
% names: node names
% Optional arguments (in order):
%   marker_size: size of node markers (default=10)

function [] = draw_graph_weight_color(vxy,edges,w,cmap_nodes, ...
    names, marker_size)

  % Optional arguments
  if (nargin < 6)
      marker_size = 10;
  end

  n = size(vxy,1);
  m = size(edges,1);
  
  vx = vxy(:,1);
  vy = vxy(:,2);

  if (isempty(w))
    w = ones(m,1);
  end

  w = w/max(abs(w));
  
  %set(gcf,'DoubleBuffer','on');
  line_width = 1;

  x1 = min(vx);
  x2 = max(vx);
  y1 = min(vy);
  y2 = max(vy);
  dx = x2-x1;
  dy = y2-y1;

  % first plot nodes
  clf;
  h=plot(vx,vy,'r.');
  axis([x1-.05*dx x2+.05*dx y1-.05*dy y2+.05*dy]);

  black = [0 0 0];
  white = [1 1 1];
  red = [1 0 0];
  blue = [0 0 1]; 
  green = [ 0 1 0];
  color1 = [1 1 1];
  color2 = [0 0 0];


  set(h,'LineWidth',line_width);
  set(h,'MarkerSize',marker_size);
  set(h,'MarkerEdgeColor',black);
  set(h,'MarkerFaceColor',black);
  axis equal;

  % then plot edges (draw weak edges first so as to not occlude stronger edges)

  ii = edges(:,1);
  jj = edges(:,2);

  absw = abs(w);
  [absw,p]=sort(absw);
  w = w(p);
  ii = ii(p);
  jj = jj(p);

  hold on;
  for k = 1:m

    i = ii(k);
    j = jj(k);

    edge = plot([vx(i) vx(j)],[vy(i) vy(j)],'-');

    set(edge,'LineWidth',w(k));
    %set(edge,'MarkerSize',marker_size);
    set(edge,'MarkerEdgeColor',black);
    set(edge,'MarkerFaceColor',black);
    if w(k) < 0
        set(edge,'Color',green);%(1+w(k))*white-w(k)*green);        
    else
        set(edge,'Color',black);%(1-w(k))*white+w(k)*blue);
    end
  end
  
  if ~isempty(cmap_nodes)
      for k = 1:n
          node = plot(vx(k), vy(k), 'ko');
          set(node, 'MarkerFaceColor', cmap_nodes(k,:));
          set(node, 'MarkerSize', marker_size);
      end
  end
      
  set(gca,'XColor','white');
  set(gca,'YColor','white');
  set(gca,'XTick',[]);
  set(gca,'YTick',[]);

  if ~isempty(names)
      for k = 1:n
          label = text(vx(k), vy(k), names{k}, 'Color', black, ...
              'VerticalAlignment','bottom', 'HorizontalAlignment','right', ...
              'FontSize', 12);
      end
  end
  
  hold off;
   
  %refresh(gcf);
  drawnow;
%  print( '-djpeg', file_name );
%  print( '-depsc', file_name );

return;


