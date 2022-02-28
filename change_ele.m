function new_arr = change_ele(arr,pos,chg)
% input:
% arr: old array
% pos: which position to change
% chg: change how much abouut this position

% output:
% new_arr is the updated new array

new_arr = arr;
new_arr(pos) = arr(pos)+chg;

end

