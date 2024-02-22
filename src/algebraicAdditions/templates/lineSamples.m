clc
clear all
close all

% Select cases
runDir = fullfile('path/to/runDir');

case1dir = string(fullfile(runDir, 'case1_dir'));
case2dir = string(fullfile(runDir, 'case2_dir'));
case3dir = string(fullfile(runDir, 'case3_dir'));
case4dir = string(fullfile(runDir, 'case4_dir'));

caseDirs = [case1dir, case2dir, case3dir, case4dir];

legends = ["Case 1", "Case 2", "Case 3", "Case 4"];

% Define lines
sameMesh = true;

L = 1;
nc = 33;
dL = L/nc;

horBegCoor = [dL/2, L/2, 0];
horEndCoor = [L-dL/2, L/2, 0];
verBegCoor = [L/2, dL/2, 0];
verEndCoor = [L/2, L-dL/2, 0];

% Plot cases
for iCase = 1:numel(caseDirs)
    caseDir = char(caseDirs(iCase));

    caseStruct = cfdReadRemoteCase(caseDir);

    if ~sameMesh || iCase==1
        horInd = cfdLineSampleIndices(caseStruct.mesh, horBegCoor, horEndCoor, nc);
        verInd = cfdLineSampleIndices(caseStruct.mesh, verBegCoor, verEndCoor, nc);
    end

    % Plot hor V
    subplot(2, 2, 1);
    hold on;
    plot(caseStruct.fluid.U.phi(horInd, 2) )

    % Plot hor p
    subplot(2, 2, 2);
    hold on;
    plot(caseStruct.fluid.p.phi(horInd))

    % Plot ver U
    subplot(2, 2, 3);
    hold on;
    plot(caseStruct.fluid.U.phi(verInd, 1))

    % Plot ver p
    subplot(2, 2, 4);
    hold on;
    plot(caseStruct.fluid.p.phi(verInd))
end

% Set titles and legends
subplot(2, 2, 1);
title('Horizontal V');
legend(legends);

subplot(2, 2, 2);
title('Horizontal p');

subplot(2, 2, 3);
title('Vertical U');

subplot(2, 2, 4);
title('Vertical p');
