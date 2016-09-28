function [MF_can, MF, FindMF] = SeekMMF(normVectors,numTrial,ConvergeAngle,ConeAngle,c,minNumSample)
% This function aims intitializing the MWO by seeking a dominat MF given a
% colleciton of surface normals via Meanshift on S2.
% output:
% MF: The non-redundant Manhattan Frame (MF) Rotation matrices.
% MF_can: MF Rotation matrices in canonical form.

MF = {};
numMF = 0;
numVectors = size(normVectors,2);
numFound = 0;
directionFound = [];
for j = 1:numTrial
    ValidMF = 0;
    
    % get random initial MF
    R = GetRandomRotation();
    
    R_update = eye(3);
    M = zeros(3,3);
    iter = 1;
%     breakF = 0;
    while acos((trace(R'*R_update) - 1)/2) > ConvergeAngle || iter == 1
        if iter ~= 1
            R = R_update;
        end
        % TODO: More smart operation should be done here. If there are only
        % two planes are observed, the other dominant direction can be
        % recovered by doing a cross product of the oether two.
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
                % breakF = 1;
                % break;
                % Here we do not straightforwardly break, if two of the
                % three direction can find the mean shift, the VO can
                % still be initialized
                if a == 2 && numFound == 0
                    break;
                end
                continue;
            end
            
            % compute mean shift
            ma = MeanShift(F,c);
            
            % compute the Ma
            alfa = norm(ma);
            ma_p = tan(alfa)/alfa*ma;
            M(:,a) = Ra' * [ma_p;1];
            M(:,a) = M(:,a) / norm(M(:,a));
            numFound = numFound + 1;
            directionFound = [directionFound a];
        end

        if numFound < 2
            numFound = 0;
            continue;
        end
        
        % if only find two dominant direction
        if numFound == 2
            v1 = M(:,directionFound(1));
            v2 = M(:,directionFound(2));
            v3 = cross(v1,v2);
            M(:,6-(directionFound(1)+directionFound(2))) = v3;
        end
        
        % force the mutual orthogonality
        [U,D,V] = svd(M);
        R_update = U*V';
        ValidMF = 1;
        iter = iter + 1;
        directionFound = [];
        numFound = 0;
    end
    
    if ValidMF == 1
        MF = [MF,{R_update}];
        numMF = numMF + 1;
        MF = MF(~cellfun('isempty',MF));
        ValidMF = 0;
    end
end

% check whether we find at least one MF
if numel(MF) == 0
    MF_can = [];
    MF = [];
    FindMF = 0;
    return;
end

% find the unique canonical form
MF_can = RemoveRedundancyMF2(MF);
FindMF = 1;
end

function R = GetRandomRotation()  
    q = rand(1,4);
    R = quat2dcm( q./norm(q) );
end

% remove the redundancy of the MF by a unique canonical form
function MF_can = RemoveRedundancyMF2(MF)%,theta)
    MF_can = {};
    % 24 possbile 
    R_poss = cell(1,24);
    cellIndex = 1;
    
    for row1=1:3
        for row2=1:3
            
            if row2 ~= row1
                for row3=1:3
                    
                    if row3 ~= row1 && row3 ~= row2
                        
                        R = zeros(3,3);
                        
                        R(1,row1) = 1.0; R(2,row2) = 1.0; R(3,row3) = 1.0;
                        if det(R) > 0 
                            R_poss{cellIndex} = R;
                            cellIndex = cellIndex + 1;
                        end
                        
                        R(1,row1) = 1.0; R(2,row2) = 1.0; R(3,row3) = -1.0;
                        if det(R) > 0 
                            R_poss{cellIndex} = R;
                            cellIndex = cellIndex + 1;
                        end
                        
                        R(1,row1) = 1.0; R(2,row2) = -1.0; R(3,row3) = 1.0;
                        if det(R) > 0 
                            R_poss{cellIndex} = R;
                            cellIndex = cellIndex + 1;
                        end
                        
                        R(1,row1) = 1.0; R(2,row2) = -1.0; R(3,row3) = -1.0;
                        if det(R) > 0 
                            R_poss{cellIndex} = R;
                            cellIndex = cellIndex + 1;
                        end
                        
                        
                        R(1,row1) = -1.0; R(2,row2) = 1.0; R(3,row3) = 1.0;
                        if det(R) > 0 
                            R_poss{cellIndex} = R;
                            cellIndex = cellIndex + 1;
                        end
                        
                        R(1,row1) = -1.0; R(2,row2) = 1.0; R(3,row3) = -1.0;
                        if det(R) > 0 
                            R_poss{cellIndex} = R;
                            cellIndex = cellIndex + 1;
                        end
                        
                        R(1,row1) = -1.0; R(2,row2) = -1.0; R(3,row3) = 1.0;
                        if det(R) > 0 
                            R_poss{cellIndex} = R;
                            cellIndex = cellIndex + 1;
                        end
                        
                        R(1,row1) = -1.0; R(2,row2) = -1.0; R(3,row3) = -1.0;
                        if det(R) > 0 
                            R_poss{cellIndex} = R;
                            cellIndex = cellIndex + 1;
                        end
                        
                    end
                    
                end
            end
        end
    end
    
    % remove redundancy and convert to canonical coordinate
    numMF = numel(MF);
    for i = 1:numMF
        minTheta = 100;
        minID = -1;
        for j = 1:24
            M_star = MF{i}*R_poss{j};
            errTheta = acos((trace(M_star)-1)/2);
            if errTheta < minTheta
                minTheta = errTheta;
                minID = j;
            end
        end
        MF_can = [MF_can,{MF{i}*R_poss{minID}}];
    end
    % clean the Null cell entries in the cell array
    MF_can = MF_can(~cellfun('isempty',MF_can));
end