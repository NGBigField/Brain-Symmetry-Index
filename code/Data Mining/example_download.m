%% Clear:
clear all; close all; clc;

%%
fileName = "myfile.mp3";
address = "www.foo.com/bar.mp3";
command_str = "curl -o "+file_name+" '"+address+"'";
% Call command:
[status,cmdout] = system(command_str,'-echo');