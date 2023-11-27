function [k_lower] = get_lower_modulation_scheme(k)
switch k
    case 6
        k_lower = 5;
    case 5
        k_lower = 4;
    case 4
        k_lower = 3;
    case 3
        k_lower = 2;
    case 2
        k_lower = 1;
    case 1
        k_lower = 0;
    otherwise
        k_lower = nan;
end
end