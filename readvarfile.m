%Function to read a .var File ,which is a SCROLL generated file.
%This function generates a 3D array image with each layer being 
%being a frame from the video .var file.
%The number of layers is the number of frames.
%
%Note: please input a filtered video file for now.
%
%Input:  flname is a filtered video .var file from SCROLL.
%
%Output: A 3D array image with each layer being 
%        being a frame from the video .var file.
%
%function [img_3D]= readvarfile(flname);
%
%e.g. img3D=readvarfile(rp02_00_flt.var);

%June 13 2012.

% modified by EAnnoni - 5/28/2015
% Adjusted X and Y to inputs for variations in video sizes
% Added lines for tmp,z,data2 to match sizes for the reshape function

function [img_3D, W, H] = readvarfile(flname)
    fid = fopen(flname, 'r', 'ieee-le.l64');  % Open the file and get the file identifier
    if fid == -1
        error('Cannot open file: %s', flname);
    end
    
    a = fread(fid);  % Read the file contents into 'a'
    fclose(fid);  % Close the file after reading

    X = a(9);  % Dimension of each frame (width)
    Y = a(5);  % Dimension of each frame (height)
    a = uint8(a);

    data_size = 1;  % Each fluorescence value is represented by 1 byte in case of a filtered file

    header_size = 24;  % Header is 24 bytes
    noofframes = (length(a) - header_size) / (X * Y) / data_size;

    data_start = header_size + 1;

    data = a(data_start:end);

    img_3D = reshape(data, [X Y ceil(noofframes)]);
    img_3D = permute(img_3D, [2 1 3]);

    W = X;
    H = Y;
end
