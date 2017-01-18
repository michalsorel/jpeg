function a =  doanalysis(z, analysis)
    if size(z,3) > 1
        no_of_frames = numel(analysis);
        a = cell(1,size(z,3));
        for i = 1:size(z,3)
            a{i} = cell(1,no_of_frames);
            for m = 1:no_of_frames
                a{i}{m} = cellfun(@(x)x/sqrt(no_of_frames),analysis{m}(z(:,:,i)),'UniformOutput', false);
            end    
        end
    else
        no_of_frames = numel(analysis);
        a = cell(1,no_of_frames);
        for m = 1:no_of_frames
            a{m} = cellfun(@(x)x/sqrt(no_of_frames),analysis{m}(z),'UniformOutput', false);
        end
    end
end