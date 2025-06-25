% Setup script for NGA-West2 correlations project
% Configures MATLAB path and environment for the reorganized project structure
%
% Run this script once after opening MATLAB to set up the workspace
% Created during project reorganization

function setup()
    % Get the current directory (should be the root project directory)
    projectRoot = pwd;
    
    % Display setup information
    fprintf('Setting up NGA-West2 correlations project...\n');
    fprintf('Project root: %s\n', projectRoot);
    
    % Define all source directories that need to be added to path
    srcDirs = {
        'src/gmpe',
        'src/residuals', 
        'src/correlations',
        'src/models',
        'src/plotting',
        'src/utils'
    };
    
    % Add each source directory to MATLAB path
    for i = 1:length(srcDirs)
        dirPath = fullfile(projectRoot, srcDirs{i});
        if exist(dirPath, 'dir')
            addpath(dirPath);
            fprintf('Added to path: %s\n', srcDirs{i});
        else
            warning('Directory not found: %s', srcDirs{i});
        end
    end
    
    % Verify key data directories exist
    dataDirs = {'data/raw', 'data/processed', 'results/figures', 'results/supplement'};
    for i = 1:length(dataDirs)
        dirPath = fullfile(projectRoot, dataDirs{i});
        if ~exist(dirPath, 'dir')
            warning('Data directory not found: %s', dataDirs{i});
        end
    end
    
    % Check for key data files
    keyDataFiles = {
        'data/raw/NGA_W2_corr_meta_data.mat',
        'data/raw/Z1_Z25_updated.mat', 
        'data/raw/durObs.mat'
    };
    
    fprintf('\nChecking for key data files:\n');
    for i = 1:length(keyDataFiles)
        filePath = fullfile(projectRoot, keyDataFiles{i});
        if exist(filePath, 'file')
            fprintf('✓ Found: %s\n', keyDataFiles{i});
        else
            fprintf('✗ Missing: %s\n', keyDataFiles{i});
        end
    end
    
    % Save the current path for future sessions (optional)
    fprintf('\nSetup complete! You can now run mainScript.m\n');
    fprintf('To save this path configuration permanently, run: savepath\n');
    
end