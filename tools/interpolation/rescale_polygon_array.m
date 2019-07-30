function result = rescale_polygon_array(x,a,b)
  concat = [];
  for i = 1:size(x,1)
    concat = [concat, x{i}];
  endfor

  %pol = x{1};
  m = min(concat');
  M = max(concat');

  for i = 1:size(x,1)
    pol = x{i};    

    if M-m<eps
        y = pol';
    else
        y = (b-a) * (pol'-m)./(M-m) + a;
    end
    
    result{i} = y';
  endfor
  
  result = result';
  
