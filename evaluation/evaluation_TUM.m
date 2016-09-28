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
load_param_MWO;

datasetName = test_data_set_Name;

% save dir
FigResultBaseDir = [inputBaseDir, datasetName, '\fig'];
if(exist(FigResultBaseDir,'dir') == 0)
    mkdir(FigResultBaseDir);
end

result_gt_gap = 3;
gt_result_gap = 0;

write2file = 1;
disp('The figures are saved in the dir:');
disp(['datasetName:',FigResultBaseDir]);

% read GT
GT = importdata('groundtruth_ass_depth_new.txt');
result = importdata([ResultBaseDir,test_data_set_Name,'\result.txt']);

numFrame = min(size(GT,1),size(result,1)); % (GT always starts after than result)
error = zeros(numFrame,1);

% rotation
roll_gt = zeros(numFrame,1);
pitch_gt = zeros(numFrame,1);
yaw_gt = zeros(numFrame,1);

roll_result = zeros(numFrame,1);
pitch_result = zeros(numFrame,1);
yaw_result = zeros(numFrame,1);

% translation
t_gt = zeros(3,numFrame);
t_result = zeros(3,numFrame);

% compute the error
t_gt0 = GT(1+gt_result_gap,2:4)';
t_result0 = result(1+result_gt_gap,2:4)';

for i = 1:numFrame
    q_gt = [GT(i+gt_result_gap,8),GT(i+gt_result_gap,5:7)];
    q_result = [result(i+result_gt_gap,8),result(i+result_gt_gap,5:7)];
    
    R_gt = quat2dcm_Eigen(q_gt);
    R_result = quat2dcm_Eigen(q_result);
    
    % load translation GT and result
    t_gt_temp = GT(i+gt_result_gap,2:4)';
    t_result_temp = result(i+result_gt_gap,2:4)';
    
    % compare rotation in terms of roll pitch yaw
    if i == 1    
        R_gt0 = R_gt;
        R_result0 = R_result;
    end
    
    error(i) = acos((trace(R_gt'*R_result)-1)/2)*57.3;
    
    [yaw_gt(i),pitch_gt(i),roll_gt(i)] = dcm2angle(R_gt);
    [yaw_result(i),pitch_result(i),roll_result(i)] = dcm2angle(R_result);
    
    if yaw_gt(i) < -0.0001
        yaw_gt(i) = yaw_gt(i) + 2*pi;
    end;
    if yaw_result(i) < -0.0001
        yaw_result(i) = yaw_result(i) + 2*pi;
    end
end

t_gt_vicon = GT(1+gt_result_gap:end,2:4)';
t_result_vicon = result(1+result_gt_gap:end,2:4)';

%% absolute error analysis
curve_x = 1:1:numFrame-3;
p = polyfit(curve_x',error(1:end-3),6);
curve_y = polyval(p,curve_x);
figure;
plot(1:1:numFrame-3,error(1:end-3),'r*');
hold on;
% plot(1:1:numFrame,error_dvo,'-b')
curve_y(curve_y < 0) = 0;
plot(curve_x,curve_y,'g-','LineWidth',5);
grid on;
title('Rotation Matrix Difference (result)');
legend('RMD','Fitting Curve','Location','northwest');
xlabel('frame number','FontSize',12,'FontWeight','bold');
ylabel('deg','FontSize',12,'FontWeight','bold');
set(gca,'FontSize',12,'FontWeight','bold');

%% Rotation
% roll
figure;
plot(1:1:numFrame,roll_gt*57.3,'-g','LineWidth',3);
hold on;
plot(1:1:numFrame,roll_result*57.3,'-r','LineWidth',3);
grid on;
title('Roll');
legend('GT','result','DVO','ICP','Location','northwest');
xlabel('frame number','FontSize',12,'FontWeight','bold');
ylabel('deg','FontSize',12,'FontWeight','bold');
set(gca,'FontSize',12,'FontWeight','bold');

% pitch
figure;
plot(1:1:numFrame,pitch_gt*57.3,'-g','LineWidth',3);
hold on;
plot(1:1:numFrame,pitch_result*57.3,'-r','LineWidth',3);
grid on;
title('Pitch');
legend('GT','result','DVO','ICP','Location','northwest');
xlabel('frame number','FontSize',12,'FontWeight','bold');
ylabel('deg','FontSize',12,'FontWeight','bold');
set(gca,'FontSize',12,'FontWeight','bold');

% yaw
figure;
plot(1:1:numFrame,yaw_gt*57.3,'-g','LineWidth',3);
hold on;
plot(1:1:numFrame,yaw_result*57.3,'-r','LineWidth',3);
grid on;
title('Yaw');
legend('GT','result','DVO','ICP','Location','northwest');
xlabel('frame number','FontSize',12,'FontWeight','bold');
ylabel('deg','FontSize',12,'FontWeight','bold');
set(gca,'FontSize',12,'FontWeight','bold');

%% translation in vicon coordinate
figure;
plot3(t_gt_vicon(1,:),t_gt_vicon(2,:),t_gt_vicon(3,:),'-g','Linewidth',3);
hold on;
plot3(t_result_vicon(1,:),t_result_vicon(2,:),t_result_vicon(3,:),'-r','Linewidth',3);
grid on;
title('Trajectory');
legend('GT','result','DVO','ICP','Location','northeast');
xlabel('X[m]','FontSize',12,'FontWeight','bold');
ylabel('Y[m]','FontSize',12,'FontWeight','bold');
zlabel('Z[m]','FontSize',12,'FontWeight','bold');
set(gca,'FontSize',12,'FontWeight','bold');
axis equal;


%% translation in x y z respectively w.r.t vicon coordinate
% translation
figure;
plot(1:1:size(t_gt_vicon,2),t_gt_vicon(1,:),'-g','LineWidth',3)
hold on;
plot(1:1:size(t_result_vicon,2),t_result_vicon(1,:),'-r','LineWidth',3);
grid on;
title('Translation in X direction');
legend('GT','result','DVO','ICP');%,'Location','northwest');
xlabel('frame number','FontSize',12,'FontWeight','bold');
ylabel('m','FontSize',12,'FontWeight','bold');
set(gca,'FontSize',12,'FontWeight','bold');
% --------------------
% set(gca,'YLim',[-1 1]);
% set(gca,'YTick',-1:0.5:1);
% ---------------------

figure;
plot(1:1:size(t_gt_vicon,2),t_gt_vicon(2,:),'-g','LineWidth',3)
hold on;
plot(1:1:size(t_result_vicon,2),t_result_vicon(2,:),'-r','LineWidth',3);
grid on;
title('Translation in Y direction');
legend('GT','result','DVO','ICP');%,'Location','northwest');
xlabel('frame number','FontSize',12,'FontWeight','bold');
ylabel('m','FontSize',12,'FontWeight','bold');
set(gca,'FontSize',12,'FontWeight','bold');

figure;
plot(1:1:size(t_gt_vicon,2),t_gt_vicon(3,:),'-g','LineWidth',3)
hold on;
plot(1:1:size(t_result_vicon,2),t_result_vicon(3,:),'-r','LineWidth',3);
grid on;
title('Translation in Z direction');
legend('GT','result','DVO','ICP','Location','northwest');
xlabel('frame number','FontSize',12,'FontWeight','bold');
ylabel('m','FontSize',12,'FontWeight','bold');
set(gca,'YLim',[0.5 2.5]);
set(gca,'YTick',0.5:0.5:2.0);
set(gca,'FontSize',12,'FontWeight','bold');