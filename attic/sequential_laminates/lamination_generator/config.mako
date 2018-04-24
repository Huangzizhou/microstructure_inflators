{
"laminates" : [
% for i in range(rank):
<%       
    scale = scales[i];
    ratio = ratios[i];
    angle = angles[i];
%>
    {
        "center" : [0, 0],
        "dimensions": [${scale}, ${ratio}],
        "rotation": ${angle}
    }
    % if i < rank-1:
    ,
    % endif
% endfor
],
"output": "rank_${rank}_laminates_${file_index}.obj",
"max_triangle_area": ${max_triangle_area},
"with_top_bottom_plates": ${with_top_bottom_plates},
"single_material": ${single_material}
}
