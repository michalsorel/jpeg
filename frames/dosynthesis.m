function s = dosynthesis(coefs, synthesis)    % assert dosynthesis(doanalysis(x,analysis), synthesis) == x
    no_of_frames = numel(synthesis);
    s = 0;
    for m = 1:no_of_frames
        s = s + synthesis{m}(coefs{m});
    end
    s = s/sqrt(no_of_frames);
end

