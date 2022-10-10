function [out] = full_nan(in)
% Convert sparse array to full, and switch zeros to nans.
out = full(in);
out(out==0) = nan;
end