function sorted_vertices = read_off_2D(filename)
% read_off - read data from OFF file.
%
%   [vertex,face] = read_off(filename);
%
%   'vertex' is a 'nb.vert x 3' array specifying the position of the vertices.
%   'face' is a 'nb.face x 3' array specifying the connectivity of the mesh.
%
%   Modified from 3D version from Gabriel Peyré
%   Copyright (c) 2003 Gabriel Peyré


fid = fopen(filename,'r');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end

str = fgets(fid);   % -1 if eof
if ~strcmp(str(1:3), 'OFF')
    error('The file is not a valid OFF one.');    
end

str = fgets(fid);
[a,str] = strtok(str); nvert = str2num(a);
[a,str] = strtok(str); nface = str2num(a);



[A,cnt] = fscanf(fid,'%f %f %f', 3*nvert);
if cnt~=3*nvert
    warning('Problem in reading vertices.');
end
A = reshape(A, 3, cnt/3);
vertex = A(1:2,:);

[A,cnt] = fscanf(fid,'%d %d %d %d\n', 3*nface);
if cnt~=3*nface
    warning('Problem in reading faces.');
end
A = reshape(A, 3, cnt/3);
face = A(2:3,:)+1;


# Sorting vertices
used_edges = false(size(face,2), 1);
polygons = {};
while sum(used_edges) != size(face,2)
  # Find unused edge
  unused_edge = -1;
  for f = 1:size(face,2)
    if !used_edges(f)
      unused_edge = f;
    endif
  endfor
  
  initial_idx = face(1,unused_edge);
  current_idx = face(2,unused_edge);
  sorted_idx = [initial_idx; current_idx];
  sorted_vertices = [vertex(:, initial_idx), vertex(:, current_idx)];
  used_edges(unused_edge) = true;

  while current_idx != initial_idx
    for f = 1:size(face,2)
      if used_edges(f)
        continue
      endif
      current_face = face(:,f);
      if current_idx == current_face(1)
        current_idx = current_face(2);
        sorted_idx = [sorted_idx; current_idx];
        sorted_vertices = [sorted_vertices, vertex(:, current_idx)];
        used_edges(f) = true;
        break
      elseif current_idx == current_face(2)
        current_idx = current_face(1);
        sorted_idx = [sorted_idx; current_idx];
        sorted_vertices = [sorted_vertices, vertex(:, current_idx)];
        used_edges(f) = true;
        break
      endif
    endfor
  endwhile  
  
  polygons = [polygons; sorted_vertices];
endwhile

sorted_vertices = polygons;

fclose(fid);
endfunction
