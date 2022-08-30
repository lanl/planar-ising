
function [] = draw_graph2(vx,vy,W)

[i,j,w] = find(W);
edges = [i,j];
draw_graph(vx,vy,edges,w);
