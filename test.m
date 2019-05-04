clc;
clear all;
close all;
disp(strcat('Corrent folder: ', pwd));


% com.mathworks.mlwidgets.html.HTMLPrefs.setUseProxy(true)
% com.mathworks.mlwidgets.html.HTMLPrefs.setProxyHost('127.0.0.1') % or whatever
% com.mathworks.mlwidgets.html.HTMLPrefs.setProxyPort('1080')


%url = 'https://lstneuro-my.sharepoint.com/personal/ariosky_areces_neuroinformatics-collaboratory_org1/_layouts/15/download.aspx?SourceUrl=%2Fpersonal%2Fariosky%5Fareces%5Fneuroinformatics%2Dcollaboratory%5Forg1%2FDocuments%2FTest%2Fdata%2Ezip';
url = 'https://drive.google.com/uc?id=1f8GCLWKbK4WpXzhESqQBFiMVNm0LulMs';
filename = strcat('data.zip');
%opts = weboptions('UserAgent','Mozilla/5.0');
options = weboptions('Timeout',Inf,'RequestMethod','get');
outfilename = websave(filename,url,options);
exampleFiles = unzip(filename,pwd);
result = true;