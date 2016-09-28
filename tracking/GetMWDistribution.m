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

% This script is used for transforming spc_MW into the 1D distribution in
% three MW direction respectively.
function D_MW = GetMWDistribution(spc_MW,SamplingInterval,scale,normalize)
    D_MW.Dx_MW = Data2OneDDistribution(spc_MW(1,:),SamplingInterval,scale,normalize);
    D_MW.Dy_MW = Data2OneDDistribution(spc_MW(2,:),SamplingInterval,scale,normalize);
    D_MW.Dz_MW = Data2OneDDistribution(spc_MW(3,:),SamplingInterval,scale,normalize);
end

% Reture 1-D density mixture distribution
function D = Data2OneDDistribution(data,SamplingInterval,scale,normalize)
    step = SamplingInterval;
    sigma = scale;
    X = min(data):step:max(data);
    numX = numel(X);
    D = zeros(2,numX);
    
    % can be faster here (lets do that in the future)
    for i = 1:numX
        D(1,i) = X(i);
        indata = data(data > X(i) - 3*sigma & data < X(i) + 3*sigma);
        D(2,i) = sum(exp(-(indata - X(i)).^2/(2*sigma^2)));
    end
    
    % Try with normalization (not sure)
    if normalize == 1
        D(2,:) = D(2,:)/max(D(2,:));
    end
end