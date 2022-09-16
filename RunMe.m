% the libraries and paths that need to be imported
% NOTE: change the paths accordingly

currPath = pwd;
basePath = '~/repos/';
libPath = [basePath 'lib/'];

% the functions
addpath(genpath([currPath '/experiments/']))
addpath([currPath '/functions/'])

% the volumes
addpath([currPath '/volumes'])

% Aspire library
cd([libPath 'ASPIRE'])
initpath;
cd(currPath)

% v1.0.1 of GlobalBioIm
% cd([libPath 'GlobalBioIm-1.0.1'])
% setMatlabPath;

% latest version of GlobalBioIm
cd([libPath 'GlobalBioIm-master'])
setGlobalBioImPath
cd(currPath)

% the LinOpCryo
addpath(genpath([currPath '/LinOpCryo']))
addpath(genpath([currPath '/LinOpCryo/functions']))
addpath(genpath([currPath '/LinOpCryo/functions_opt']))

% create folder to save intermediate results
%if ~isfolder('./interm_res') % recent matlab version
if ~exist('./interm_res', 'dir')
	mkdir('./interm_res/')
end
