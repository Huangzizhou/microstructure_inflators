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
"output": "rank_2_laminates_${triangle_density}.obj",
"max_area": ${max_area}
}
