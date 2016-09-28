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
clc;
clear;
close all;
% load parameter and dir
addpath('..\');
addpath('..\evaluation');
load_param_MWO;

% Outer loop (data)
inputDataDir = [inputBaseDir, test_data_set_Name];
ResultBaseDir = [ResultBaseDir, test_data_set_Name];
SpcDir = [ResultBaseDir,'\SPC\'];
if( exist(inputDataDir,'dir') == 0 )
    mkdir(ResultBaseDir);
end

% load input
% read the image (timestamp) name list
depthMapFiles = dir([inputDataDir,'\depth']);
depthMapList = depthMapFiles(3:end); % ignore the first two entries which are . and ..
numFrames = size(depthMapList,1);

% output file
if saveResult == 1
    fileID = fopen([ResultBaseDir,'\result.txt'],'w');
end

% groundtruth of the first frame (for converting the result to w.r.t vicon)
q_gt0 = [0.2239,-0.4871,0.7673,-0.3519];
R_gt0 = quat2dcm_Eigen(q_gt0);
t_gt0 = [-2.5508;0.9872;1.1019];

% main loop
for i = 1:numFrames   
    % directly load the sn and spc
    data = importdata([SpcDir,depthMapList(i).name,'.mat']);
    sn = data.sn;
    spc = data.PointCloud;
    
    % we only use points whose depth is between max(0.5,d_min)
    % and 2*median(d) - d_min
    d_min = max(min(spc(3,:)),0.5);
    d_med = median(spc(3,:));
    spc = spc(:,spc(3,:) > d_min & spc(3,:) < 2*d_med - d_min);
    
    % Tracking
    % Motion notation: R is from Manhattan Frame (MF) to camera coordinate
    %                  t is from MF to camera coordinate
    if InitMfFound == 1 % see whether the MF has been found
        % Tracking MF
        [R,IsTracked] = TrackingMF(R,sn,ConvergeAngle,ConeAngle_tracking,c,minNumSample);

        % if lost tracking
        if IsTracked == 0
            disp(num2str(i));
            disp('lost tracking!');
            break;
        end
        
        % compensate the rotation by convert the point cloud to MF frame.
        spc_MW2 = R'*spc;
        
        % translation estimation
        D_MW_old = D_MW1;
        [t_r, D_MW1] = EstimateTranslation(D_MW1,spc_MW2,config,options);
        D_MW_new = D_MW1;
        t = t - t_r;
    else
        % Initialization (Seek the dominant MF)
        [MF_can,MF,FindMF] = SeekMMF(sn,numTrial,ConvergeAngle,ConeAngle,c,minNumSample);
        if(FindMF == 1)
            R = ClusterMMF(MF_can,ratio);% The ouput is a cell of several MF_nonRd
            if isempty(R) == 0
                InitMfFound = 1;
                disp('Initialization done!');
            end
            R = R{1};

            % compensate the rotation by convert the point cloud to MF frame.
            spc_MW1 = R'*spc;
            D_MW1 = GetMWDistribution(spc_MW1,config.SamplingInterval,config.scale,config.normalize);           
            
            % Initial translation is zero
            t = [0 0 0]';% the translation (position) is defined in MW frame
            
            % for result convertion
            if convertVicon == 1
                R_result0 = R;
                t_result0 = t;
            end
        end
    end
    
    % record pose estimation
    if InitMfFound == 1       
        % make the result consistent with groundtruth, which are the
        % rotation from camera to vicon (R_v_c) and the location of camera w.r.t
        % vicon (t_v_c). Users need to manually copy the groundtruth of the
        % first frame to here, otherwise the system will output the result
        % w.r.t the first frame.
        if convertVicon == 1
            % rotation
            R_wrt_vicon = R_gt0*R_result0*R';
            q_wrt_vicon = dcm2quat_Eigen(R_wrt_vicon);
            %translation
            t_result_delta = t - t_result0;
            t_wrt_vicon = t_gt0 + R_gt0*R_result0*t_result_delta;
        else
            % convert R to quartnion q
            q = dcm2quat_Eigen(R);
        end
        
        % record the estimated pose
        if saveResult == 1
            if convertVicon == 0
                fprintf(fileID,'%s %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f\n',num2str(depthMapList(i).name(1:17)),t(1),t(2),t(3),q(2),q(3),q(4),q(1));
            else
                fprintf(fileID,'%s %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f\n',num2str(depthMapList(i).name(1:17)),t_wrt_vicon(1),t_wrt_vicon(2),t_wrt_vicon(3),q_wrt_vicon(2),q_wrt_vicon(3),q_wrt_vicon(4),q_wrt_vicon(1));
            end
        end
    end
    disp(['The ',num2str(i),' frame']);
end

if saveResult == 1
    fclose(fileID);
end

%% evaluate
evaluation_TUM;