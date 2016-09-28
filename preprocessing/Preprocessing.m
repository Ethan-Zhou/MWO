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

% This is script aims at creating the super point cloud (3D patch) in advance for efficient experiment.
clc;
clear;
close all;

% load parameter and dir
addpath('..\');
load_param_MWO;

% load input
inputDataDir = [inputBaseDir,test_data_set_Name];
ResultDir = [ResultBaseDir,test_data_set_Name];
if(exist(ResultDir,'dir'))
    mkdir(ResultDir);
end

disp(['Processing: ',inputDataDir]);
depthMapFiles = dir([inputDataDir,'\depth']);
depthMapNameList = depthMapFiles(3:end);% the 1st and 2nd element are '.' and '..'
numFrames = size(depthMapNameList,1);
saveDir = [ResultDir,'\SPC\'];
if(exist(saveDir,'dir') == 0)
    mkdir(saveDir);
end

% main loop
for i = 1:numFrames
    % load the depthMap
    depthMap = imread([inputDataDir,'\depth\',depthMapNameList(i).name]);
    
    % surface normal fitting
    [sn,spc,PointCloud] = GetSurfaceNormalCell(depthMap,cellsize,K);
    
    % save the spc and sn
    save([saveDir,depthMapNameList(i).name,'.mat'],'sn','spc','PointCloud');
    
    % display processing status
    disp(num2str(i));
end