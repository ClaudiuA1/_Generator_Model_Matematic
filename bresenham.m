%% Funcția Bresenham (pentru a trasa o linie între două puncte de pe grilă)
function coords = bresenham(x1, y1, x2, y2)
% bresenham - Returnează coordonatele celulelor de pe linia dintre (x1,y1) și (x2,y2)
dx = abs(x2 - x1);
dy = abs(y2 - y1);
sx = sign(x2 - x1);
sy = sign(y2 - y1);
err = dx - dy;

x = x1;
y = y1;
coords = [x, y];
while ~(x == x2 && y == y2)
    e2 = 2 * err;
    if e2 > -dy
        err = err - dy;
        x = x + sx;
    end
    if e2 < dx
        err = err + dx;
        y = y + sy;
    end
    coords(end+1,:) = [x, y];  %#ok<AGROW>
end
end
