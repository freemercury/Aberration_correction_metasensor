% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Preprocessing and reaglin for raw scanning lightfield data.
%  The Code is created based on the method described in the following paper 
%   [1]  JIAMIN WU, YUDUO GUO and CHAO DENG etc,
%        An integrated imaging sensor for aberration-corrected 3D photography
%        Nature, 2022. 
% 
%   Contact: YUDUO GUO (gyd20@mails.tsinghua.edu.cn)
%   Date  : 6/10/2022
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfolder_code = pwd;
cd('..');
cfolder = pwd;
para.resize = 0.9908;
para.rotate = 0;
para.raw_folder = [cfolder,'/Data/Rawdata'];
para.save_folder = [cfolder,'/Data/Rawdata/Calibrate_Raw/'];
mkdir(para.save_folder)
cd(cfolder_code);
Output_name = Meta_Calibrate(para.raw_folder,[1],para.resize,para.rotate,[para.save_folder]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
autoCenterMode = 0;                             % Auto Find the center point
centerX=3894;                                  
centerY=2979;  
startFrame = 0;                                 % Start Frame: 0 ~ N-1
frameJump  = 1;                                 % Frame jump step
centerView  = [cfolder,'/Data/Rawdata/Meta_image/centerView.tiff'];
                                                % '.None.tiff','.only.tiff'
realignName = [cfolder,'/Data/Rawdata/Meta_image/test'];% Realign data name
mkdir([cfolder,'/Data/Rawdata/Meta_image/']);
Nnum = 15;
Nshift = 5;                                     % Shift number, [1,3,5,15]
Nx = 250;                                       % Half Microlens number picked along X
Ny = 190;                                        % Half Microlens number picked along Y
groupCount = 1;                                 % Group number of Meta image
groupMode = 0;                                  % 0 = (1-9,10-18,...);1 = (1-9,2-10,...) 
autoCenterFrame = startFrame;                                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name = Output_name;                             % The first Tiff data name, normally '**.0.tiff'    
confName = [cfolder,'/Code/Realign/',num2str(Nshift),'x',num2str(Nshift),'.conf.sk.png'];
                                                % Scanning configuration
resize = 0;                                     % Resize to 15X15 if need
rotation =  0;                                  
preResize = 1.0;                                
slightRotation = 0;                            
warpFile = 'warp.tiff';                         
flatFile = 'flat.tiff';
if(autoCenterMode == 1)
    img = double(imread(name, autoCenterFrame + 1));
    [centerX, centerY, ~] = AutoCenter(img, Nnum);
end
cd ./Realign
command = sprintf('ReAlign %d %s %s %d %d %d %d %d %d %s %d %d %d %d %s %s %d %f %f %s %s',...
Nshift, realignName, name, startFrame, centerX, centerY, Nx,...
groupCount, groupMode, confName, frameJump, resize, Nnum, Ny, 'ZGX', centerView, rotation, preResize, slightRotation, warpFile, flatFile);
system(command);

function [Xcenter,Ycenter, prob] = AutoCenter(img, Nnum)

    SUB_NUM = 3;%n
    SUB_MAX_RANGE_RATIO = 5;%k
    MAX_AREA_TRUST_RATIO = 4;%m
    img = double(img);
    fullX = size(img, 2);
    fullY = size(img, 1);
    %整体，最亮处(kxk pix)(m倍可信度)，分块nxn
    ansList = zeros(1 + MAX_AREA_TRUST_RATIO + SUB_NUM*SUB_NUM, 3);
    fprintf('Full Image :');
    [ansList(1,1), ansList(1,2), ansList(1,3)] = AutoCenterSub(img, Nnum, 0, 0);
    fprintf('\n');
    maxImg = max(img(:));
    [maxPosY, maxPosX] = find(img == maxImg);
    maxPosX = maxPosX(round(size(maxPosX, 1) / 2));
    maxPosY = maxPosY(round(size(maxPosY, 1) / 2));
    
    rangeSX = round(maxPosX - fullX / SUB_MAX_RANGE_RATIO);
    if(rangeSX < 1) rangeSX = 1; end
    rangeEX = round(maxPosX + fullX / SUB_MAX_RANGE_RATIO);
    if(rangeEX > fullX) rangeEX = fullX; end
    
    rangeSY = round(maxPosY - fullY / SUB_MAX_RANGE_RATIO);
    if(rangeSY < 1) rangeSY = 1; end
    rangeEY = round(maxPosY + fullY / SUB_MAX_RANGE_RATIO);
    if(rangeEY > fullY) rangeEY = fullY; end
    fprintf('Lightest x%d:', MAX_AREA_TRUST_RATIO);
    [ansList(2,1), ansList(2,2), ansList(2,3)] = ...
        AutoCenterSub(img(rangeSY : rangeEY, rangeSX : rangeEX), Nnum, rangeSX - 1, rangeSY - 1);
    fprintf('\n');
    ansList(2:2+MAX_AREA_TRUST_RATIO-1, :) = repmat(ansList(2,:), [MAX_AREA_TRUST_RATIO, 1]);
    
    anchorPtX = zeros(SUB_NUM+1, 1);
    for i = 0:SUB_NUM
        anchorX = round(1 + i * fullX/SUB_NUM);
        if(anchorX > fullX) anchorX = fullX; end
        anchorPtX(i+1) = anchorX;
    end
    anchorPtY = zeros(SUB_NUM+1, 1);
    for i = 0:SUB_NUM
        anchorY = round(1 + i * fullY/SUB_NUM);
        if(anchorY > fullY) anchorY = fullY; end
        anchorPtY(i+1) = anchorY;
    end
    
    idx = 1 + MAX_AREA_TRUST_RATIO + 1;
    
    for i = 1:SUB_NUM%x
        for j = 1:SUB_NUM%y
            fprintf('Sub X=%d Y=%d:', i,j);
            [ansList(idx,1), ansList(idx,2), ansList(idx,3)] = ...
                AutoCenterSub(img(anchorPtY(j) : anchorPtY(j+1), anchorPtX(i) : anchorPtX(i+1)), Nnum, anchorPtX(i) - 1, anchorPtY(j) - 1);
            distance = sqrt((i - ((SUB_NUM + 1) /2))^2 + (j - ((SUB_NUM + 1) /2))^2);
            prob_loss = 1 - (distance / SUB_NUM);
            fprintf('prob loss ratio = %.2f\n',  prob_loss);
            ansList(idx,3) = ansList(idx,3) * prob_loss;
            idx = idx + 1;
        end
    end
    savedAnsList = ansList;
    ansList(:,3) = ansList(:,3) / sum(ansList(:,3));
    ansList(:,1) = ansList(:,1) .* ansList(:,3);
    ansList(:,2) = ansList(:,2) .* ansList(:,3);
    myAns = round([sum(ansList(:,1)), sum(ansList(:,2))]);
    prob = 0;
    probCount = 0;
    for i = 1:size(savedAnsList, 1)
        if(myAns == savedAnsList(i, 1:2))
            probCount = probCount + 1;
            prob = prob + savedAnsList(i, 3);
        end
    end
    if(probCount ~= 0)
        prob = prob / probCount;
    end
    Xcenter = myAns(1) + floor(size(img,2) / Nnum / 2) * Nnum;
    Ycenter = myAns(2) + floor(size(img,1) / Nnum / 2) * Nnum;
    fprintf('AutoCenter found x = %d, y = %d (range from 0~[size-1]), credibility = %.2f%%, check at ImageJ\n', Xcenter, Ycenter, prob*100);
end

function [Xcenter, Ycenter, prob] = AutoCenterSub(img, Nnum, x_offset, y_offset)
    
    BEST_RATIO = 0.3;
    WORST_RATIO = 0.9;

    img = double(img);
    img = img ./ max(img(:));
    img = img.^2;
    kernal = fspecial('gaussian',[Nnum,Nnum],3);
    img = imfilter(img, kernal);
    locMatrix = zeros(Nnum, Nnum);
    for i = 1:Nnum
        for j = 1:Nnum
            picMat = img(i:Nnum:end, j:Nnum:end);
            %locMatrix(i,j) = mean(mean(picMat));
            avg = mean(mean(picMat));
            picMat(picMat < avg) = 0;
            hugePos = find(picMat ~= 0);
            locMatrix(i,j) = sum(picMat(:)) / size(hugePos, 1);
        end
    end
    
    sumX = sum(locMatrix);
    sumY = sum(locMatrix,2);
    sumX = sumX + circshift(sumX, 1, 2);
    sumY = sumY + circshift(sumY, 1, 1);
    darkX = 1;
    for i = 1:Nnum
        if(sumX(i) < sumX(darkX))
            darkX = i;
        end
    end
    darkY = 1;
    for i = 1:Nnum
        if(sumY(i) < sumY(darkY))
            darkY = i;
        end
    end
    Xcenter = mod(darkX + floor(Nnum / 2) - 1 + x_offset, Nnum);
    Ycenter = mod(darkY + floor(Nnum / 2) - 1 + y_offset, Nnum);
    prob = (WORST_RATIO - min(min(locMatrix)) / max(max(locMatrix))) / (WORST_RATIO - BEST_RATIO);
    if(prob > 1) 
        prob = 1; 
    elseif(prob < 0) 
        prob = 0; 
    end
    fprintf('AutoCenterSub x = %d, y = %d, prob = %f ', Xcenter, Ycenter, prob);
end