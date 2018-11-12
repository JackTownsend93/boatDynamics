% Test for xml read/write.

xDoc    = xmlread('params.xml');
params  = xmlwrite(xDoc);

params = strsplit(params, '\n')';

% Change speed (inlet velocity) through iterations.
changeSpeed = input('Would you like to change SPEED over time? (y/n):\n', 's');
% Change geometry through iterations.
changeGeom  = input('Would you like to change GEOMETRY over time? (y/n):\n', 's');


%if changeSpeed == 'y'
%    speedInit  = input('What is the initial speed (m/s)?\n');
%    speedFinal = input('What is the final speed (m/s)?\n');


% SPEED.

% If changing speed.
if changeSpeed == 'y'
    speed = input('What speed?\n');
    % Write speed and reference speed into params.
    speedString = strcat('      <inletVelocity>', num2str(speed,'%.3f'), '</inletVelocity>');
    refSpeedString = strcat('      <uRef>', num2str(speed,'%.3f'), '</uRef>');
    params(25,1) = {speedString};
    params(26,1) = {refSpeedString};
    % Save to xml file.
    fid = fopen('params.xml','wt');
    for i = 1:length(params)
        fprintf(fid, '%s\n', params{i});
    end    
    fclose(fid);
    
    % Notify.
    fprintf('Speed is variable.\n');

% If not changing speed.
else if changeSpeed == 'n'
        % Notify.
        fprintf('Speed is invariable.\n');
    
    % Handle input errors.
    else
    fprintf('You must enter y/n.\n');
    end
end

% GEOMETRY.





