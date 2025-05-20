function [letter_list,ip_letter_list] = ans_w08_readdata(input_letter)
%% Function ans_w08_readdata
%   for mcb111 homework 08
%   take in a character, for example, 'A' or 'B'
%   return two matrices, reference for perfect letter and imperfect letter
%% Perfect reference
filename = strcat('abc.',input_letter,'.cat');
str = fileread(filename); % read perfect reference as a long string
parts = strtrim(regexp( str, '(\n\n)+', 'split')); % split into letters
letter_list = zeros(25, length(parts)); % matrix to store all letters
for i = 1:length(parts) % start processing each letter
    lines = strtrim(regexp(  parts{i}, '(\r|\n)+', 'split')); %split lines in each letter
    letter = (cell2mat(lines) =='x')*2 - 1; %convert into -1,1 vector
    letter_list(:, i) = letter; % save
end
%% imperfect refrences:
ip_filename = strcat('abc.pmin0.3.pmax0.5.',input_letter,'.dat');
ip_str = fileread(ip_filename); % read imperfect as a long string
ip_parts = reshape(ip_str,30,[]); % splt into letters
ip_letter_list = zeros(25, length(ip_parts)); % matrix to store all imperfect letters
for i =1:length(ip_parts) % start processing each imperfect letter!
    ip_lines = strtrim(regexp(ip_parts(:,i)', '(\r|\n)+', 'split')); %split lines in each letter
    ip_letter =(cell2mat(ip_lines)=='x')*2 - 1; %convert into -1,1 vector
    ip_letter_list(:,i) = ip_letter; % save
end
end

