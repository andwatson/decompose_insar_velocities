function [comp_perp1,comp_perp2] = null_line(comp1,comp2)
%=================================================================
% function null_line()
%-----------------------------------------------------------------
% Calculates the null line from two line-of-sights (i.e. the line along
% which two LOS are completely insensitive to ground motion.
% Returns the two planes for decomposition.
%
%                                                                  
% INPUT:                                                           
%   comp1: [E N U] component vectors for first look direction
%   comp2: [E N U] component vectors for second look direction
% OUTPUT:    
%   comp_perp1: plane 1 (perpendicular to null line)
%   comp_perp2: plane 2
%   
% Andrew Watson     13-06-2022
%                                                                  
%=================================================================

% find component vectors of null line
comp_null = cross(comp1,comp2);

% back calculate inc and az out of interest
inc_null = acosd(comp_null(3));
az_null = acosd(comp_null(1)./sind(inc_null));

% find orthogonal vectors
orth = null(comp_null(:).');
comp_perp1 = orth(:,1);
comp_perp2 = orth(:,2);

end