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

function [t, D_MW2] = EstimateTranslation(D_MW1,spc_MW2,config,options)
% This version is based on discrete distribution correlation (22/12/2015).
% input: 
%       spc_MW: is the super point cloud described in the MW frame, which
%       mean the rotation has be compensated.
%       
% convert spc to three 1D discrete distributions (make it parallel)
Dx_MW1 = D_MW1.Dx_MW;
Dy_MW1 = D_MW1.Dy_MW;
Dz_MW1 = D_MW1.Dz_MW;
[Dx_MW2,Dy_MW2,Dz_MW2] = Spc2OneDDistribution(spc_MW2,config.SamplingInterval,config.scale,config.normalize);
D_MW2.Dx_MW = Dx_MW2;
D_MW2.Dy_MW = Dy_MW2;
D_MW2.Dz_MW = Dz_MW2;

% do 1D align for each channel (make it parallel in the future)
config_x = initialize_config(Dx_MW1,Dx_MW2,config);
if config.DisplayDistribution == 1
    DisplayEnergy(config_x);
end
sx = DistributionAlign(config_x,options);

config_y = initialize_config(Dy_MW1,Dy_MW2,config);
if config.DisplayDistribution == 1
    DisplayEnergy(config_y);
end
sy = DistributionAlign(config_y,options);

config_z = initialize_config(Dz_MW1,Dz_MW2,config);
if config.DisplayDistribution == 1
    DisplayEnergy(config_z);
end
sz = DistributionAlign(config_z,options);

t = [sx;sy;sz];
end

% Return 2*M matrix, the first dimension stores the sampling position, the
% second dimension stores the intensity
function [Dx, Dy, Dz] = Spc2OneDDistribution(spc_MW,SamplingInterval,scale,normalize)
    Dx = Data2OneDDistribution(spc_MW(1,:),SamplingInterval,scale,normalize);
    Dy = Data2OneDDistribution(spc_MW(2,:),SamplingInterval,scale,normalize);
    Dz = Data2OneDDistribution(spc_MW(3,:),SamplingInterval,scale,normalize);
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
    
    % Try with normalization
    if normalize == 1
        D(2,:) = D(2,:)/max(D(2,:));
    end
end

% Nonlinear optimization
function t = DistributionAlign(config,options)
   t = fmincon(@costFunc, config.init_t, [ ], [ ], [ ], [ ], config.Lb, config.Ub, [ ], options, config);
end

% CorrelationDistance
function f = CorrelationDistance(t,config)
    A_t = config.A;
    A_t(1,:) = A_t(1,:) + t;
    B = config.B;
    SamplingInterval = config.SamplingInterval;

    % smarter
    % 1. find the entries in A_t (the transformed one) that are locating within the range of B
    idA_inB = find(A_t(1,:) > B(1,1) & A_t(1,:) < B(1,end));
    idA = 1:1:size(A_t,2);
    idA_outB = setdiff(idA,idA_inB);
    
    % 2. the intensity of the ones locating outside the range of B will be
    % directly add to the cost
    if config.distMeasure == 2
        f = sum(A_t(2,idA_outB)*A_t(2,idA_outB)');%lets try SAD measure
    end
    if config.distMeasure == 1
        f = sum(abs(A_t(2,idA_outB)));
    end
    
    % 3. create the interpolated vector, and calculate the other part of
    % the cost
    A_t_inB = A_t(:,idA_inB);
    length_A_t_inB = size(A_t_inB,2);
    n = floor((A_t_inB(1,1)-B(1,1))/SamplingInterval);
    id = n+1;
    B1 = B(:,id:id+length_A_t_inB-1);
    B2 = B(:,id+1:id+length_A_t_inB);
    d2 = B(1,1)+(n+1)*SamplingInterval-A_t(1,1);
    alpha = d2/SamplingInterval;
    B_p = alpha*B1(2,:) + (1-alpha)*B2(2,:);
    % get distance
    if config.distMeasure == 2
        f = sum((A_t_inB(2,:)-B_p)*(A_t_inB(2,:)-B_p)') + f;
    end
    % try SAD
    if config.distMeasure == 1
        f = sum(abs(A_t_inB(2,:)-B_p)) + f;
    end
end

% costFunc
function [f,g] = costFunc(t,config)
    f = CorrelationDistance(t,config);
    f1 = CorrelationDistance(t-config.dt,config);
    f2 = CorrelationDistance(t+config.dt,config);
    g = (f2 - f1)/(2*config.dt);
end

% 
function configExtend = initialize_config(model, scene, config)
configExtend = config;
configExtend.B = scene;
maxBx = max(scene(1,:));
minBx = min(scene(1,:));
model = model(:,model(1,:) > minBx & model(1,:) < maxBx);
lengthA = size(model,2);
truncatedLength = floor(lengthA*0.05);% truncate A
configExtend.A = model(:,truncatedLength+1:lengthA-truncatedLength);
end

%
function DisplayEnergy(configObj)
cost_x = -0.1:0.001:0.1;
cost = zeros(1,numel(cost_x));
for i = 1:numel(cost_x)
    cost(i) = costFunc(cost_x(i),configObj);
end
figure;
plot(cost_x,cost,'g-','LineWidth',3);
hold on;
[mini,miniId] = min(cost);
plot(cost_x(miniId),cost(miniId),'ro','LineWidth',3);
hold off;

figure;
subplot(1,2,1);
stem(configObj.A(1,:),configObj.A(2,:),'*');
hold on;
stem(configObj.B(1,:),configObj.B(2,:),'o');

subplot(1,2,2);
stem(configObj.A(1,:)+cost_x(miniId),configObj.A(2,:),'*');
hold on;
stem(configObj.B(1,:),configObj.B(2,:),'o');
hold off;

figure;
stem(configObj.A(1,1:5:end),configObj.A(2,1:5:end),'*');
hold on;
stem(configObj.B(1,1:5:end),configObj.B(2,1:5:end),'o');
axis off;
% export_fig distribution_alignment.png -transparent
end
