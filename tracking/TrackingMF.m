%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:   Yi Zhou                                                          *
% Contact:  yi.zhou@anu.com                                                  *
% License:  Copyright (c) 2016 Yi Zhou, ANU. All rights reserved.            *
%                                                                            *
% Redistribution and use the code, with or without                           *
% modification, are permitted provided that the following conditions         *
% are met:                                                                   *
% * Redistributions of source code must retain the above copyright           *
%   notice, this list of conditions and the following disclaimer.            *
% * Neither the name of ANU nor the names of its contributors may be         *
%   used to endorse or promote products derived from this software without   *
%   specific prior written permission.                                       *
%                                                                            *
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"*
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE  *
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE *
% ARE DISCLAIMED. IN NO EVENT SHALL ANU OR THE CONTRIBUTORS BE LIABLE        *
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL *
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT         *
% LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY  *
% OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF     *
% SUCH DAMAGE.                                                               *
%*****************************************************************************/

function [R_update,updated,denMedian] = TrackingMF(R,normVectors,ConvergeAngle,ConeAngle,c,minNumSample)
% This function aims at tracking the MF.
% Output:
%       R_update: the current orientation of the sensor

R_update = R';%infact, R_updata can be randomly assigned
updated = 0;
numFound = 0;
directionFound = [];
% uncertain threshold
iter = 1;
% density temp
denTemp = zeros(3,1) + 0.00001;
while acos((trace(R'*R_update)-1)/2) > ConvergeAngle || iter == 1
        if iter ~= 1
            R = R_update;
        end
        for a = 1:3
            Ra = [R(:,mod(a+3,3)+1),R(:,mod(a+4,3)+1),R(:,mod(a+5,3)+1)]';
            F = [];
            % Let's make it fast!
            nVps = Ra*normVectors;
            eta = sqrt(nVps(1,:).*nVps(1,:) + nVps(2,:).*nVps(2,:));
            index = find(eta < sin(ConeAngle));
            nVps_inlier = nVps(:,index);
            tan_alfa = eta(index)./abs(nVps(3,index));
            alfa = asin(eta(index));
            f = [alfa./tan_alfa.*nVps_inlier(1,:)./nVps_inlier(3,:);
                 alfa./tan_alfa.*nVps_inlier(2,:)./nVps_inlier(3,:)];
             
            select = ~isnan(f);
            select2 = select(1,:).*select(2,:);
            select3 = find(select2 == 1);
            F = f(:,select3);
            
            if numel(F) < minNumSample
                % Here we do not straightforwardly break, if two of the
                % three direction can find the mean shift, the VO can
                % still be initialized
                if a == 2 && numFound == 0
                    % lost tracking
                    R_update = [];
                    updated = 0;
                    denMedian = 0;
                    % visualize the surface normal
%                     vx = normVectors(1,:);
%                     vy = normVectors(2,:);
%                     vz = normVectors(3,:);
%                     figure;
%                     plot3(vx,vy,vz,'g.');
%                     axis off;
%                     axis equal;
                    return;
                end
                continue;
            end
            
            % compute mean shift
            [ma,denTemp(a,1)] = MeanShift(F,c);
            
            % compute the Ma
            alfa = norm(ma);
            ma_p = tan(alfa)/alfa*ma;
            R_update(:,a) = Ra' * [ma_p;1];
            R_update(:,a) = R_update(:,a) / norm(R_update(:,a));
            numFound = numFound + 1;
            directionFound = [directionFound a];
        end
        
        if numFound < 2
            % lost tracking
            R_update = [];
            updated = 0;
            denMedian = 0;
%             % visualize the surface normal
%                     vx = normVectors(1,:);
%                     vy = normVectors(2,:);
%                     vz = normVectors(3,:);
%                     figure;
%                     plot3(vx,vy,vz,'g.');
%                     axis off;
%                     axis equal;
%                     return;
            return;
        end
        
        % if only find two dominant direction
        if numFound == 2
            v1 = R_update(:,directionFound(1));
            v2 = R_update(:,directionFound(2));
            v3 = cross(v1,v2);
            R_update(:,6-(directionFound(1)+directionFound(2))) = v3;
        end
        
        % force the mutual orthogonality (each direction with different weight based on the local density)
        denMedian = median(denTemp);% for debug
%         denTemp(min_den_id,1) = 1;% let the least un-trustable direction with weight 1
        [U,D,V] = svd([R_update(:,1)*denTemp(1,1),R_update(:,2)*denTemp(2,1),R_update(:,3)*denTemp(3,1)]);
        R_update = U*V';
        iter = iter + 1;
        directionFound = [];
        numFound = 0;
end
updated = 1;