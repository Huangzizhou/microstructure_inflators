function im = off2png(mesh_path, image_path)
  V = read_off_2D(mesh_path);
  V = rescale_polygon_array(V, .00, 1.0);
  
  N = size(V, 1);
  
  % Fix all polygons that don't have correct orientation
  for i = 1:N
    polygon = V{i};
    
    vx = polygon(1, :);
    vy = polygon(2, :);

    x = [1.1];
    y = [1.1];

    in = inpolygon(x, y, vx, vy);

    if sum(in) == 1
        V{i} = fliplr(polygon)
    endif
  endfor

  % Find boundary
  boundary_idx = 0;
  for i = 1:N
    candidate = V{i};
    
    vx = candidate(1, :);
    vy = candidate(2, :);

    x = [];
    y = [];
    for j = 1:N
        if i == j
            continue
        endif

        test_polygon = V{j};
        x = [x; test_polygon(1,1)];
        y = [y; test_polygon(2,1)];
    endfor

    in = inpolygon(x, y, vx, vy);

    if sum(in) == (N-1)
        boundary_idx = i;
        break
    endif
    %pause
  endfor


  t = linspace(0.0, 1.0, 201);
  [x,y] = meshgrid(t,t);
  
  result = ones(201,1);  
  for i = 1:N
    polygon = V{i};
    vx = polygon(1, :);
    vy = polygon(2, :);
    
    in = inpolygon (x, y, vx, vy);

    if i != boundary_idx
        in = ~in;
    endif

    result = result & in;
    
    %plot (vx, vy, x(in), y(in), "r+", x(!in), y(!in), "bo");
    %axis ([0, 1, 0, 1]); 
    %pause;
  endfor
    
  %plot (vx, vy, x(result), y(result), "r+", x(!result), y(!result), "bo");
  %axis ([0, 1, 0, 1]); 


  im = reshape(result, 201, 201);
  imwrite(~im, image_path);
  
endfunction
