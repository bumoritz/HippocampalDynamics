function [varargout]=paq2lab_ov20220105(varargin)
% function [varargout]=paq2lab(varargin);
% Read binary files acquired with PackIO into MATLAB workspace
% Return recording of each channel as a vector array
% Based on paqread (used in DaqViewer)
% Usage: [data, names, units, rate,filename, fposition]=paq2lab;

% Adam Packer
% January 30th, 2007

%MB20210615: commented out from line 112

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose and open file

if ~isempty(varargin);%if a file path specified, use it
    fullpath = varargin{1};
    [pathstr, name, ext] = fileparts(fullpath);
    filename = name;
else%if no path already, ask the user for it
    [filename,pathname,FilterIndex]=uigetfile({'*.bin;*.paq'},'Choose a data file');
    if ~FilterIndex
        varargout(1)={[]};
        varargout(2)={[]};
        varargout(3)={[]};
        varargout(4)={[]};
        varargout(5)={[]};
        varargout(6)={[]};
        return
    end
    fullpath=[pathname filename];
end

fid=fopen(fullpath);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in rate, number of channels, and channel names
rate=fread(fid,1,'float32','b');
numchans=fread(fid,1,'float32','b');
for i=1:numchans;
    number_of_characters=fread(fid,1,'float32','b');
    channelname{i}=[];
    for j=1:number_of_characters
        channelname{i}=[channelname{i}, strrep(fread(fid,1,'float32=>char', 'b'),' ','')];
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in hardware channel ('HWchan') and units if available (*.paq only)
[pathstr, name, ext] = fileparts(fullpath);
if strcmp(ext,'.paq')
    for k=1:numchans
        number_of_characters=fread(fid,1,'float32','b');
        HWchan{k}=[];
        for m=1:number_of_characters
            HWchan{k}=[HWchan{k}, strrep(strrep(fread(fid,1,'float32=>char','b'),' ',''),'/','_')];
        end
    end
    for n=1:numchans
        number_of_characters=fread(fid,1,'float32','b');
        units{n}=[];
        for q=1:number_of_characters
            units{n}=[units{n}, strrep(fread(fid,1,'float32=>char','b'),' ','')];
        end
    end
elseif strcmp(ext,'.bin')
    for k=1:numchans
        HWchan{k}='unknown';
        units{k}='unknown';
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get filesize and fposition (which should be after all header info)
dirinfo=dir(pathstr);
for n=1:size(dirinfo);
    if strcmp(dirinfo(n).name,[name ext]);
        filesize=dirinfo(n).bytes;
    end
end

fposition=ftell(fid);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get 'info' or 'channels' (default)
try
    secondinput=lower(varargin{2});
catch
    secondinput='channels';
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch secondinput
    case 'channels'
        % Try to set ephysdefaults, or ask no questions
        try
            [rangequest,qmquestion,chunkquest,channels,filter,stackquest,laserquest]=ephysdefaults;
        catch
            rangequest=0;
            qmquestion=0;
            channels=1;
        end
        
        % Try to get included channels, or ask user
        try
            includedchannels=varargin{3};
        catch
            includedchannels=1:length(channelname); %MB20210615
%             if channels==1;
%                 for a=1:numchans
%                     if strcmp(ext,'.paq');
%                         charcell{a}=strcat(HWchan{a},'-',channelname{a});
%                     end
%                     if strcmp(ext,'.bin');
%                         charcell{a}=strcat('0','-',channelname{a});
%                     end
%                 end
% 
%                 [includedchannels,ok]=listdlg('ListString',charcell,'PromptString','Plot which channels... All?',...
%                     'Name','Channel Selection','InitialValue',1:length(charcell));%find out which channels to plot
%                 if isempty(ok) || ~ok;
%                     varargout(1)={[]};
%                     varargout(2)={[]};
%                     varargout(3)={[]};
%                     varargout(4)={[]};
%                     varargout(5)={[]};
%                     varargout{6}={[]};
%                     return
%                 end
%             elseif channels==0;
%                 includedchannels=1:length(channelname);
%             end
        end          
        if rangequest==1;
            prompt1=['Filesize(MB) = ',...
                sprintf('%-.2f\n',filesize/10^6),...
                '1sec data (MB) = ',...
                sprintf('%-.2f\n',4*rate*numchans/10^6),...
                'Filesize(sec) = ',...
                sprintf('%-.2f\n',(filesize/10^6)/(4*rate*numchans/10^6)),...
                'Read from time(sec)'];
            prompt={prompt1 'Stop reading at time(sec) [Enter zero for all]'};
            dlg_title=['Select range to import:'];
            num_lines=1;
            def={'0','0'};
            answer=APinputdlg(prompt,dlg_title,num_lines,def);
            if isempty(answer)
                varargout(1)={[]};
                varargout(2)={[]};
                varargout(3)={[]};
                varargout(4)={[]};
                varargout(5)={[]};
                varargout{6}={[]};
                return
            end
            starttime=str2num(cell2mat(answer(1)));
            stoptime=str2num(cell2mat(answer(2)));
            startsample=starttime*rate;
            stopsample=stoptime*rate;
            startbyte=startsample*numchans;
            stopbyte=stopsample*numchans;
            if (stopbyte > filesize) ||(stopbyte==0)
                stopbyte=filesize;
            end

            if qmquestion==0;
                button='Quick';
            elseif qmquestion==1;
                button = questdlg('Read quickly OR Open more samples?','Reading style','Quick','More','Quick');
            end
            if isempty(button)
                varargout(1)={[]};
                varargout(2)={[]};
                varargout(3)={[]};
                varargout(4)={[]};
                varargout(5)={[]};
                varargout{6}={[]};
                return
            elseif strcmp(button,'More')        % READ BY FOR LOOP THROUGH DATA
                waithandle = waitbar(0,'Reading from PackIO file');%for user
                for i=1:length(includedchannels)
                    fseek(fid,fposition+startbyte*4+includedchannels(i)*4-4,'bof');
                    data(:,i)=fread(fid,(stopbyte-startbyte)/numchans,'*float32',4*(numchans-1),'b');
                    waithandle = waitbar(i/length(includedchannels),waithandle);
                end
                close(waithandle);
                fclose(fid);
            elseif strcmp(button,'Quick')   % READ BY AUTO COLUMN SORT FREAD
                fseek(fid,fposition+startbyte*4,'bof');
                data=fread(fid,[numchans,(stopbyte-startbyte)/numchans],'*float32','b');
                notchans=1:numchans;
                notchans(includedchannels)=[];
                data(notchans,:)=[];
                data=data';
                fclose(fid);
            end
        elseif rangequest==0;
            if qmquestion==0;
                button='Quick';
            elseif qmquestion==1;
                button = questdlg('Read quickly OR Open more samples?','Reading style','Quick','More','Quick');
            end
            if isempty(button)
                varargout(1)={[]};
                varargout(2)={[]};
                varargout(3)={[]};
                varargout(4)={[]};
                varargout(5)={[]};
                varargout{6}={[]};
                return
            elseif strcmp(button,'More')        % READ BY FOR LOOP THROUGH DATA
                waithandle = waitbar(0,'Reading from PackIO file');%for user
                for i=1:length(includedchannels)
                    fseek(fid,fposition+includedchannels(i)*4-4,'bof');
                    data(:,i)=fread(fid,filesize/numchans,'*float32',4*(numchans-1),'b');
                    waithandle = waitbar(i/length(includedchannels),waithandle);
                end
                close(waithandle);
                fclose(fid);
            elseif strcmp(button,'Quick')   % READ BY AUTO COLUMN SORT FREAD
                fseek(fid,fposition,'bof');
                data=fread(fid,[numchans,filesize/numchans],'*float32','b');
                notchans=1:numchans;
                notchans(includedchannels)=[];
                data(notchans,:)=[];
                data=data';
                fclose(fid);
            end
        end
    case 'info'
        % Read in the Object information and Event information.
        info.ObjInfo.SampleRate=rate;
        info.HwInfo='Probably acquired with PackIO';
        info.ObjInfo.SamplesAcquired=((filesize-fposition)/4)/numchans;
        for a=1:numchans;
            if strcmp(ext,'.paq')
                info.ObjInfo.Channel(a).HwChannel=HWchan{a};
                info.ObjInfo.Channel(a).Units=units{a};
            else
                info.ObjInfo.Channel(a).HwChannel=0;
                info.ObjInfo.Channel(a).Units='SeePackIO';
            end
            info.ObjInfo.Channel(a).ChannelName=channelname{a};
        end
        data = info;
        fclose(fid);
        includedchannels=[];
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
varargout(1)={data};
varargout(2)={channelname(includedchannels)};
varargout(3)={units(includedchannels)};
varargout(4)={rate};
varargout(5)={filename};
varargout{6}={fposition};