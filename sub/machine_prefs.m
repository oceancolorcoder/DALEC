function [fontName, projectPath, imageryPath] = machine_prefs()
% function [fontName, projectPath, imageryPath] = machine_prefs()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismac
    fontName = 'Times New Roman';
    projectPath = '/Users/daurin/Projects/';
    imageryPath = '/Volumes/Megalodon';
elseif isunix
    fontName = 'Times';
    projectPath = '/accounts/daurin/Projects/';
    imageryPath = '/glusteruser/daurin/';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%