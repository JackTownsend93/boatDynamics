function speedPrevious = funcParamsChange(changeSpeed, speedPrevious, speedIncrement)

% This function, if the user has asked to vary speed, writes a new speed to the 
% params.xml file on restarting each sim.
% The speed written in this iteration is passed on in the main script loop to be
% used to calculate the next value (speedPrevious + speedIncrement).

% NOTES: - Strictly reads from and writes to the file named "params.xml" in the working directory.
%        - Can not handle any comment lines in the params.xml file.
%        - Assumes certain values (e.g.: inletVelocity) are on specific lines in the params.xml
%          file. Ensure the file matches these assumptions.



if changeSpeed == 'y'
    % Calc new speed.
    speed = speedPrevious + speedIncrement;
    
    % Read in the params.xml file as cell array, delimiting by returns:
    xDoc   = xmlread('params.xml');
    params = xmlwrite(xDoc);
    params = strsplit(params, '\n')';

    % Write to appropriate cells (inletVel and uRef).
    speedString = strcat('      <inletVelocity>', num2str(speed,'%.3f'), '</inletVelocity>');
    refSpeedString = strcat('      <uRef>', num2str(speed,'%.3f'), '</uRef>');
    params(25,1) = {speedString};
    params(26,1) = {refSpeedString};

    % Save as an xml file in current directory.
    fid = fopen('params.xml','wt');
    for i = 1:length(params)
        fprintf(fid, '%s\n', params{i});
    end    
    fclose(fid);

    % Save speed to be used as speedPrevious for next iteration.
    speedPrevious = speed;

    % Notify.
    fprintf('%s: Speed incremented for next sim.\n', mfilename);

elseif changeSpeed == 'n'
    fprintf('%s: No speed change.\n', mfilename); 

end

end