function cells = floodFill(mat, i, j, rows, cols)
stack = [i, j]; cells = []; visited = false(rows, cols);
while ~isempty(stack)
    [ci, cj] = deal(stack(end,1), stack(end,2));
    stack(end,:) = [];
    if visited(ci,cj), continue; end
    visited(ci,cj) = true;
    if ~strcmp(mat{ci,cj}, '1'), continue; end
    cells(end+1,:) = [ci, cj];
    for d = [-1 0; 1 0; 0 -1; 0 1]'
        ni=ci+d(1); nj=cj+d(2);
        if ni>=1 && ni<=rows && nj>=1 && nj<=cols && ~visited(ni,nj)
            stack(end+1,:) = [ni, nj];
        end
    end
end
end