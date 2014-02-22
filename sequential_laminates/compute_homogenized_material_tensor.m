%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute_homogenized_material_tensor.m
% 02/20/2014 - Qingnan Zhou
% (partially extracted from sequentialLaminates.m by Julian)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the homogenized material property with the specified parameters.
% This method is partially extracted and refectored from sequentialLaminates.m
% for debugging and prototyping purposes.  It is NOT optimized for speed or
% memory.  This code works in 2D only!
%
% @param[in] lamA, muA  Lame parameters for material A
% @param[in] lamB, muB  Lame parameters for material B
% @param[in] p          number of lamination steps
% @param[in] ratios     An array of ratios of material A in each lamination
%                       layer
% @param[in] angles     An array of anlges in degrees.  Each angle specifies the
%                       angle between lamination direction with the Y axis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A_star = compute_homogenized_material_tensor(lamA, muA, lamB, muB, p, ratios, angles)
    assert(p == size(ratios, 2));
    assert(p == size(angles, 2));
    A = [lamA+2*muA, lamA,       0;
         lamA,       lamA+2*muA, 0;
         0,          0,          2*muA];
    B = [lamB+2*muB, lamB,       0;
         lamB,       lamB+2*muB, 0;
         0,          0,          2*muB];
    function fA = fAGen(e)
        e1e1 = e(1) * e(1); e2e2 = e(2) * e(2); e1e2 = e(1) * e(2);
        fA = (1 / (4 * muA)) * [4*e1e1, 0,           2*e1e2;
                                0,           4*e2e2, 2*e1e2;
                                2*e1e2,      2*e1e2, e1e1+e2e2] + ...
            (1 / (2 * muA + lamA) - 1 / muA) * [e1e1*e1e1, e1e1*e2e2, e1e1*e1e2;
                                                e1e1*e2e2, e2e2*e2e2, e1e2*e2e2;
                                                e1e1*e1e2, e1e2*e2e2, e1e1*e2e2];
    end

    e_array = ones(2,p);
    for eIt = 1:p
        e_array(1,eIt) = cos(degtorad(angles(eIt)));
        e_array(2,eIt) = sin(degtorad(angles(eIt)));
    end

    ratio_prod = 1.0;
    for itr = 1:p
        ratio_prod = ratio_prod * (1.0 - ratios(itr));
    end

    tmp_matrix = (B - A)^(-1);
    for itr = 1:p
        tmp_ratio_prod = ratios(itr);
        for jtr =1:itr-1
            tmp_ratio_prod = tmp_ratio_prod * (1.0 - ratios(jtr));
        end
        tmp_matrix = tmp_matrix + tmp_ratio_prod * fAGen(e_array(:,itr));
    end

    A_star = A + ratio_prod * tmp_matrix^(-1);

end
