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

% InvPbLib
cd([libPath 'InvPbLib'])
setMatlabPath;
cd(currPath)

% the LinOpCryo
addpath(genpath([currPath '/LinOpCryo']))
addpath(genpath([currPath '/LinOpCryo/functions']))
addpath(genpath([currPath '/LinOpCryo/functions_opt']))

% create folder to save intermediate results
if ~isfolder('./interm_res')
	mkdir('./interm_res/')
end
