%Title : DvDt max and Activation time for cells
%Created by: Sanaz Hosseini
%Created:November 2023

%% Set up
clc, clear all, close all;
% Read the .var file
[mov, W, H] = readvarfile('pacing02_1000S_filter.var'); %M4TRI1_rep2_20231004
mov = double(mov);
%Selects a region of interest (ROI) from the video frames (80x80 pixels)
videoWidthnew = 80; %Extracts color (RGB) and grayscale intensity information from the frames.
videoHeightnew = 80;
numFrames = 4997;
totalTime = 10; % seconds
frameRate = 500; % frames per second
mov=mov(1:videoWidthnew,1:videoHeightnew,:); 
 for i=1:1:numFrames              
                   shot=mov(:,:,i);
                    
                    imagesc(shot)
          colormap(gray(256))
          colorbar
                    caxis([0 200]);
                    F(i) = getframe(gcf)
                   
 end
% Create a VideoWriter object to save the video
  writerObj = VideoWriter('pacing02_1000S_filter.avi');
  writerObj.FrameRate = 500; %stem cell
  %writerObj.FrameRate = 1000; % real heart
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);

% Initialize arrays to store color and intensity information
colorData = zeros(videoHeightnew, videoWidthnew, numFrames, 3); % RGB color
intensityData = zeros(videoHeightnew, videoWidthnew, numFrames); % Grayscale intensity

% Read and process each frame from the .var file
for frameNumber = 1:numFrames
    frame = mov(:, :, frameNumber); % Read the frame from the .var file

    % Store color information (RGB)
    colorData(:, :, frameNumber, 1) = frame;

    % Convert the frame to grayscale and store intensity information
    intensityData(:, :, frameNumber) = frame;
end

%% Filter the noise with the Mask

% Calculate the time interval (in milliseconds) between frames
mSeconds = 1000 / frameRate;
% Create a directory for saving figures
%outputDirectory = '/Users/Sanaz/Desktop/video/figs';
%mkdir(outputDirectory); % Create the directory if it doesn't exist

MaskThresholdEditFieldValue = 15;
Mask = zeros(80, 80) + 1;
maximum = max(max(max(mov)));

Thresh = (MaskThresholdEditFieldValue / 100) * maximum;

for i = 1:80
    for j = 1:80
        pixelIntensity = squeeze(intensityData(i, j, :));

        if max(pixelIntensity) >= Thresh
            Mask(i, j) == 1
            % Display the figure and perform additional actions
            numFramesThisPixel = length(pixelIntensity); % Get the actual number of frames for this pixel
            timeVector = (1:numFramesThisPixel) * mSeconds;
            % figure;
            % First subplot for the Mask filter
            % subplot(2, 1, 1);
            %  plot(timeVector, pixelIntensity);
            %  xlabel('Time (ms)');
            %  ylabel('Intensity');
            %  title(['Pixel Intensity at (' num2str(i) ',' num2str(j) ') - Mask Filter']);

            % Apply moving average filter
            windowSize = 5;
            filteredIntensity = movmean(pixelIntensity, windowSize);

            % Second subplot for the moving average filter
            %  subplot(2, 1, 2);
            %  plot(timeVector, filteredIntensity);
            %  xlabel('Time (ms)');
            %  ylabel('Filtered Intensity');
            %  title(['Filtered Pixel Intensity at (' num2str(i) ',' num2str(j) ') - Moving Average Filter']);
            % % Save the figure as an image file with both subplots
            %  filename = fullfile(outputDirectory, ['second filter _Pixel_' num2str(i) '_' num2str(j) '.png']);
            %  saveas(gcf, filename);
            %  close(gcf);
        else
          Mask(i, j) = 0;
        end
    end
end


%% Calculate max dv/dt
startTime = 1900;  % Start time in milliseconds
endTime = 2500;    % End time in milliseconds

% Convert the time interval to frame indices
startFrame = round(startTime / mSeconds);
endFrame = round(endTime / mSeconds);

% Initialize arrays to store the maximum dv/dt and the time for each pixel
maxDvDt = zeros(80, 80);
timeOfMaxDvDt = zeros(80, 80);

% Determine the total number of pixels in the ROI
totalPixels = sum(Mask(:));

% Initialize the matrix to store pixel information
pixelInfoMatrix = zeros(totalPixels, 3); % Assuming 80x80 pixels

% Initialize pixel index
pixelIndex = 1;

% Loop through each pixel in the ROI
for i = 1:80
    for j = 1:80
        if Mask(i, j) == 1 
            % Extract pixel intensity data for the chosen pixel
            pixelIntensity = squeeze(intensityData(i, j, startFrame:endFrame));

            % Apply moving average filter
            windowSize = 5;
            filteredIntensity = movmean(pixelIntensity, windowSize);
        
            % Calculate the derivative of intensity over time (dv/dt)
            dvDt = diff(filteredIntensity) / mSeconds;
        
            % Find the maximum dv/dt and the time at which it occurs
            [maxDvDt(i, j), maxIndex] = max(dvDt);
        
            % Calculate the time at which the maximum dv/dt occurs
            timeOfMaxDvDt(i, j) = startTime + maxIndex * mSeconds;

            % Store information in the matrix
            pixelInfoMatrix(pixelIndex, 1) = sub2ind([80, 80], i, j); % Pixel number
            pixelInfoMatrix(pixelIndex, 2) = maxDvDt(i, j);             % Max dv/dt
            pixelInfoMatrix(pixelIndex, 3) = timeOfMaxDvDt(i, j);       % Time of Max dv/dt

            % Increment the pixel index
            pixelIndex = pixelIndex + 1;

            % Print the results for each pixel
            fprintf('Pixel (%d, %d): Max dv/dt = %.2f, Time = %.2f ms\n', i, j, maxDvDt(i, j), timeOfMaxDvDt(i, j));
        end
    end
end


%%
% Create a colormap plot for max dv/dt
figure;
h = imagesc(maxDvDt);

% Set white background for the colormap where dv/dt is zero
colormap('parula'); % can change 'parula' to any other colormap you prefer
set(gca, 'Color', 'w'); % Set the background color to white

% Replace zero values with NaN to distinguish from actual data
maxDvDt(maxDvDt == 0 ) = NaN;
maxDvDt(maxDvDt == 0) = NaN;

maxDvDt(maxDvDt <=  0.100000000000001 | maxDvDt>=10) = NaN;


%Remove max dv/dt for pixels x=1 and x=5
%maxDvDt (53:58 , 74:80)= NaN
%maxDvDt (80 , :) = NaN;
%maxDvDt(36:40, 77:80) = NaN;
%maxDvDt(1:5, 76:80) = NaN;
%maxDvDt ( 80, 1:16) = NaN
%maxDvDt (58, 23)  = NaN
%maxDvDt (59, 23) = NaN
%maxDvDt(49,74) = NaN;
%maxDvDt(17:18, 79:80) = NaN;
% Save the figure as an image file (e.g., PNG)
saveas(gcf, 'max_dvdt_plot.png');

colorbar;
title('Maximum dv/dt');
xlabel('Pixel X');
ylabel('Pixel Y');

% Correct the Y-axis orientation
axis xy;

% Set white background for regions where dv/dt is zero,
set(h, 'AlphaData', ~isnan(maxDvDt));

% Reverse the dv/dt plot from top to bottom
set(gca, 'YDir', 'reverse');

% Adjust color limits for better visualization
caxis([0, max(maxDvDt(:))]);

% Add more ticks to the colorbar with a range
ticks = 0:0.3:max(maxDvDt(:));
colorbar('Ticks', ticks);

% Manually set y-axis tick locations at intervals of 10 from 0 to the max value
yticklocations = 0:10:80;
yticks(yticklocations);

% Reverse the y-axis tick labels
yticklabels_reversed = flip(cellstr(num2str(yticklocations')));
set(gca, 'YTickLabel', yticklabels_reversed); % Set the reversed labels

% Save the figure as an image file (e.g., PNG)
saveas(gca, 'max_dvdt_plot.png');

%%
% Create a scatter plot for activation time
figure;
% Replace zero values with NaN to distinguish from actual data
%timeOfMaxDvDt(timeOfMaxDvDt == 0 | timeOfMaxDvDt >2600) = NaN;
timeOfMaxDvDt(timeOfMaxDvDt== 0) = NaN;
%timeOfMaxDvDt ( 53:58, 74:80) = NaN;
%timeOfMaxDvDt (1:5, 76:80) = NaN;
%timeOfMaxDvDt ( 80, 1:16) = NaN
% timeOfMaxDvDt  (49,74) = NaN;
% timeOfMaxDvDt (17:18, 79:80) = NaN;
timeOfMaxDvDt (1, 1) = NaN;
% Create meshgrid for x and y
[x, y] = meshgrid(1:size(maxDvDt, 2), 1:size(maxDvDt, 1));

% Scatter plot without reversing x-axis
scatter(x(:), y(:), 50, timeOfMaxDvDt(:), 'filled');
title('Activation Time (AT)');
xlabel('Pixel X');
ylabel('Pixel Y');
colorbar;

% Set color bar limits for the scatter plot with an upper limit of 2500
startTime = 1900; 
caxis([startTime, 2500]); % upper limit as needed

colormap('parula'); 
% Set white background for the scatter plot where either dv/dt or time is zero
set(gca, 'Color', 'w'); % Set the background color to white

% Manually set x-axis tick locations and labels at intervals of 10 from 0 to 80
xticklocations = 0:10:80;
xticklabels = xticklocations;
xticks(xticklocations);

% Add more ticks to the colorbar with a range of 100 between each two ticks
ticks = startTime:50:2500; 
colorbar('Ticks', ticks);

% Reverse the dv/dt plot from top to bottom
set(gca, 'YDir', 'reverse');

% Manually set y-axis tick locations and labels at intervals of 10 from 0 to the max value
yticklocations = 0:10:80;
yticks(yticklocations);

% Reverse the y-axis tick labels
yticklabels_reversed = flip(cellstr(num2str(yticklocations')));
set(gca, 'YTickLabel', yticklabels_reversed); % Set the reversed labels

% Save the figure as an image file (e.g., PNG)   
saveas(gcf, 'activation_time_plot.png');
