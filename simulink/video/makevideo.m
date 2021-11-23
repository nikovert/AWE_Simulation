function makevideo(output_name, preamble)

    if nargin <2
        preamble = 'video/vidPic_';
    end
    if nargin < 1
        output_name = 'pumping_cycle';
    end
       
    writerObj = VideoWriter(output_name,'MPEG-4');

    open(writerObj);
    
    K = 1;
    datenum_prev = 0;
    
    while K <= inf
        filename = sprintf([preamble,'%d.png'], K);
        if ~isfile(filename)
            break
        end
        FileInfo = dir(filename);
        if FileInfo.datenum < datenum_prev 
            break
        else
            datenum_prev = FileInfo.datenum;
        end
        thisimage = imread(filename);
        writeVideo(writerObj, thisimage);
        K = K +1; 
    end
    close(writerObj);
end