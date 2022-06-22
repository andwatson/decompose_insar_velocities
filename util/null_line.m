function [comp_perp1,comp_perp2] = null_line(comp1,comp2)
%=================================================================
% function null_line()
%-----------------------------------------------------------------
% Calculates the null line from two line-of-sights.
%                                                                  
% INPUT:                                                           
%   los1: [E N U]
%   los2: [E N U]
% OUTPUT:    
%   track_vel: 
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