% Jason K. Johnson
% LANL, Aug 2010
%
% SYNOPSIS draw an undirected graph with edge weights controlling
% shading of drawn edges.
% INPUTS:
% vx,vy: horizontal and vertical locations of vertices in the plane
% edges: m x 2 matrix of edges {edges(k,1),edges(k,2)}, k = 1,...,m
% w: edge weights

function [] = draw_graph(vxy,edges,w,file_name)

  n = size(vxy,1);
  m = size(edges,1);
  
  vx = vxy(:,1);
  vy = vxy(:,2);

  if (isempty(w))
    w = -ones(m,1);
  end
  
  w(w>0) = w(w>0)/max(w);
  %w = w/max(w);

  %set(gcf,'DoubleBuffer','on');
  line_width = 2;
  marker_size = 16;

  x1 = min(vx);
  x2 = max(vx);
  y1 = min(vy);
  y2 = max(vy);
  dx = x2-x1;
  dy = y2-y1;

  % first plot nodes
  clf;
  h=plot(vx,vy,'ko');
  label_posx = vx + [0 0 .05*dx 0 .03*dx]';
  label_posy = vy + [.1*dy .1*dy .05*dy .1*dy .05*dy]';
  text(label_posx, label_posy, {'a', 'd', 'c', 'b', 'e'}, 'FontSize', 24)
  axis([x1-.05*dx x2+.05*dx y1-.05*dy y2+.05*dy]);

  black = [0 0 0];
  white = [1 1 1];
  blue = [0 0.2 0.5]; 


  set(h,'LineWidth',line_width);
  set(h,'MarkerSize',marker_size);
  set(h,'MarkerEdgeColor',black);
  set(h,'MarkerFaceColor',black);
  axis equal;

  % then plot edges (draw weak edges first so as to not occlude stronger edges)

  ii = edges(:,1);
  jj = edges(:,2);

  [w,p]=sort(w);
  ii = ii(p);
  jj = jj(p);

  hold on;
  for k = 1:m
    if w(k) >= -.5

        i = ii(k);
        j = jj(k);

        edge = plot([vx(i) vx(j)],[vy(i) vy(j)],'ko-');

        set(edge,'LineWidth',line_width);
        set(edge,'MarkerSize',marker_size);
        set(edge,'MarkerEdgeColor',black);
        set(edge,'MarkerFaceColor',black);
        if w(k) < 0
            w(k) = 0;
        elseif w(k) > 1
            w(k) = 1;
        end
        weight1 = exp(w(k)) / (exp(w(k)) + exp(1-w(k)));
        weight2 = exp(weight1) / (exp(weight1) + exp(1-weight1));
        weight = power(w(k),1/1.7);
        set(edge,'Color',(1-weight)*white+weight*black);
%        set(edge,'Color',(1-w(k))*white+w(k)*black);
    end
  end
  
  for k = 1:m

    if w(k)< -.5
        i = ii(k);
        j = jj(k);

        edge = plot([vx(i) vx(j)],[vy(i) vy(j)],'ko-');

        set(edge,'LineWidth',line_width);
        set(edge,'MarkerSize',marker_size);
        set(edge,'MarkerEdgeColor',black);
        set(edge,'MarkerFaceColor',black);
%        set(edge,'Color','red');
    end
  end


  set(gca,'XColor','white');
  set(gca,'YColor','white');
  set(gca,'XTick',[]);
  set(gca,'YTick',[]);

  hold off;
   
  %refresh(gcf);
  drawnow;
%  print( '-djpeg', file_name );
%  print( '-depsc', file_name );
print(file_name, '-dpdf')

return;


