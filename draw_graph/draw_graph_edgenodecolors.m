% Jason K. Johnson
% LANL, Aug 2010
%
% SYNOPSIS draw an undirected graph with edge weights controlling
% shading of drawn edges.
% INPUTS:
% vx,vy: horizontal and vertical locations of vertices in the plane
% edges: m x 2 matrix of edges {edges(k,1),edges(k,2)}, k = 1,...,m
% w: edge weights

function [] = draw_graph_edgenodecolors(vxy,edges,w,w_nodes,file_name, names)

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
  marker_size = 6;

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

    edge = plot([vx(i) vx(j)],[vy(i) vy(j)],'r.-');

    set(edge,'LineWidth',line_width);
    set(edge,'MarkerSize',marker_size);
    set(edge,'MarkerEdgeColor',black);
    set(edge,'MarkerFaceColor',black);
    if w(k) < 0
        set(edge,'Color',green);%(1+w(k))*white-w(k)*green);        
    else
        set(edge,'Color',black);%(1-w(k))*white+w(k)*blue);
    end
  end
  
  if ~isempty(w_nodes)
%       w_nodes = w_nodes/max(abs(w_nodes));

      for k = 1:n
          %if(strcmp(char(w_nodes(k)),'D'))%
              if ( 0.5 < w_nodes(k) && w_nodes(k) < 1.5 ) %Democrats
            node=plot(vx(k),vy(k),'bo');
            set(node,'MarkerFaceColor',blue);%(1+w_nodes(k))*color1);
          %elseif(strcmp(char(w_nodes(k)),'R'))%
              elseif ( 1.5 < w_nodes(k) && w_nodes(k) < 2.5 ) %Republicans
            node=plot(vx(k),vy(k),'ro');
            set(node,'MarkerFaceColor',red);%(1+w_nodes(k))*color1);
          else %Others
            node=plot(vx(k),vy(k),'ko');
            set(node,'MarkerFaceColor',black);%(1-w_nodes(k))*color2);
          end
      end
  end
      
  set(gca,'XColor','white');
  set(gca,'YColor','white');
  set(gca,'XTick',[]);
  set(gca,'YTick',[]);

  if ~isempty(names)
      for k = 1:n
          if(strcmp(char(w_nodes(k)), 'D')) % Democrats
              label = text(vx(k), vy(k), names{k}, 'Color', blue, ...
                  'VerticalAlignment','bottom', 'HorizontalAlignment','right');
          elseif(strcmp(char(w_nodes(k)), 'R')) % Republicans
              label = text(vx(k), vy(k), names{k}, 'Color', red, ...
                  'VerticalAlignment','bottom', 'HorizontalAlignment','right');
          else % Others
              label = text(vx(k), vy(k), names{k}, 'Color', black, ...
                  'VerticalAlignment','bottom', 'HorizontalAlignment','right');
          end
      end
  end
  
  hold off;
   
  %refresh(gcf);
  drawnow;
%  print( '-djpeg', file_name );
%  print( '-depsc', file_name );

return;


