clc
clear all
close all
fps=500;

[mov,W,H] = readvarfile('M4EPDFB_1_0807_pacing02_1000Sfilter.var'); %%load data 
mov = double(mov);


    %% Masking

                MaskThresholdEditFieldValue=20;
                        Mask = zeros(80,80)+1;
                        maximum =  max(max(max((mov))))
                        Thresh = (MaskThresholdEditFieldValue/100) * maximum; % the highest masked value
                        for i = 1:80
                            for j = 1:80
                                m = squeeze(mov(i,j,:));
                                if max(m) < Thresh
                                    Mask(i,j) = 0;
                                end
                            end
                        end                       
%                         Mask(1:20,:)=0; %%Modify if needed

                    
figure
                       imagesc(Mask)

 %% APD80 Map
 StartFrame=1287; %EPC_1 0720  Choose one peak
EndFrame=2784;

 mov = mov(:,:,StartFrame:EndFrame);
            L = length(mov);
            APD_map = zeros(W,H);
            iteration = 1;
            APDpercent = 80;
            APDpercent = APDpercent/100;
            window_left = 250;  
            window_right = 250; 
            numberpeaks = 1;
            
            
            for i = 1 :1: 80
                i
                for j = 1 :1:80
                    
ensemble_pixel_data=[];


                    g = squeeze(mov(i,j,:));
                    B = 1/10*ones(10,1);
                    g = filter(B,1,g);
                    plot(g)
                    thres_h = APDpercent*max(g); % Set thresh_h to APD% of the max for pixel (i,j)
                    [pks,locs] = findpeaks(g,'MinPeakHeight',thres_h,'MinPeakProminence',5); %find local maxima above thresh_h
                   
                    if (isempty(pks) | thres_h < 20 | length(pks) < numberpeaks) %If no peaks were found or thresh_h was small or only only two or few peaks were found
                        APD_map(i,j) = NaN;
                        continue;
                    end
                    count = 1;
                    for k = 1 : length(locs)
                        l_limit = locs(k) - window_left;
                        u_limit = locs(k) + window_right;
                     if (l_limit > 1 & u_limit < L) % if 110 frames to the left and right of the peak is still part of the movie:
                            ensemble_pixel_data(count,:)=g(l_limit:u_limit); %set the count(th) row of enseble_pixel_data to be the frame range 110 frames before and after a peak for pixel (i,j)
                            count = count+1;
                        end
                    end

                    if isempty(ensemble_pixel_data)
                        APD_map(i,j) = NaN;
                        continue;
                    end


                    if size(ensemble_pixel_data,1)>1 %if the number of rows in the ensemble pixel data matrix (number of pixels with peaks)
                        ensemble_pixel_data = mean(ensemble_pixel_data);
                    end
                    APD_80_threshold = (max(ensemble_pixel_data)-(APDpercent*(max(ensemble_pixel_data) - min(ensemble_pixel_data))));
                    


                    APD_80_low = max(find(ensemble_pixel_data(1:window_left+5) < APD_80_threshold));
                    if isempty(APD_80_low)
                        APD_80_low = 1;
                    end
                    APD_80_high = min(find(ensemble_pixel_data(window_left+1:end) < APD_80_threshold)) + window_left;
                    if isempty(APD_80_high)
                        APD_80_high = window_left+window_right;
                    end
                    APD_80 = APD_80_high - APD_80_low;
                    APD_80_ms = APD_80/fps;
                    APD_map(i,j) = APD_80_ms.*1000; %convert to milliseconds
                    ensemble_pixel_data = [];


                end
            end
            
APD_map=APD_map.*Mask;

imagesc(APD_map)
colorbar
caxis([0 1000])
        myColorMap = parula(1000);                    
            myColorMap(1,:) = 1;           
            colormap(myColorMap);
           title('APD80 Color Map')