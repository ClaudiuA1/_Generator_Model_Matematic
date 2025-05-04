% --- getAdjacentNodes helper ---
function nodes = getAdjacentNodes(nodeMap,i,j,rows,cols)
nodes=[];
if i>1, nodes(end+1)=nodeMap(i-1,j); end
if i<rows, nodes(end+1)=nodeMap(i+1,j); end
if j>1, nodes(end+1)=nodeMap(i,j-1); end
if j<cols, nodes(end+1)=nodeMap(i,j+1); end
nodes = unique(nodes(nodes>0));
end