function varargout = lineage_launcher_REALTIME_fish_3(varargin)
%lineage_launcher_REALTIME_fish_3 M-file for lineage_launcher_REALTIME_fish_3.fig
%      lineage_launcher_REALTIME_fish_3, by itself, creates a new lineage_launcher_REALTIME_fish_3 or raises the existing
%      singleton*.
%
%      H = lineage_launcher_REALTIME_fish_3 returns the handle to a new lineage_launcher_REALTIME_fish_3 or the handle to
%      the existing singleton*.
%
%      lineage_launcher_REALTIME_fish_3('Property','Value',...) creates a new lineage_launcher_REALTIME_fish_3 using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to lineage_launcher_REALTIME_fish_3_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      lineage_launcher_REALTIME_fish_3('CALLBACK') and lineage_launcher_REALTIME_fish_3('CALLBACK',hObject,...) call the
%      local function named CALLBACK in lineage_launcher_REALTIME_fish_3.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to runbutton (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lineage_launcher_REALTIME_fish_3

% Last Modified by GUIDE v2.5 06-Oct-2016 10:07:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lineage_launcher_REALTIME_fish_3_OpeningFcn, ...
                   'gui_OutputFcn',  @lineage_launcher_REALTIME_fish_3_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before lineage_launcher_REALTIME_fish_3 is made visible.
function lineage_launcher_REALTIME_fish_3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no outputname args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line outputname for lineage_launcher_REALTIME_fish_3
handles.outputname = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes lineage_launcher_REALTIME_fish_3 wait for user response (see UIRESUME)
% uiwait(handles.figure1);
handles.ROIcount=0;
handles.ROIs={};
guidata(hObject,handles)


% --- Outputs from this function are returned to the command line.
function varargout = lineage_launcher_REALTIME_fish_3_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning outputname args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line outputname from handles structure
varargout{1} = handles.outputname;


% --- Executes on button press in runbutton.
function runbutton_Callback(hObject, eventdata, handles)
% hObject    handle to runbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
suffix=get(handles.embryosuffix,'String');
parameterfilenames={};

embryodir=handles.embryodir;
wrongslashes=findstr(embryodir,'\');
embryodir(wrongslashes)='/';
outputdir=get(handles.outputdir,'String');
if(iscell(outputdir))
    outputdir=outputdir{1};
end
wrongslashes=findstr(outputdir,'\');
outputdir(wrongslashes)='/';
outputdir=[outputdir,'/'];


outputdirectory=[outputdir,handles.embryodirname,'/'];
mkdir(outputdirectory);

if length(handles.ROIs)>1
    handles.ROIs = {handles.ROIs{1}};
end

for i=1:length(handles.ROIs)
%for online analysis discard all but the first ROI
%for each embryo ROI
%make copy of param file with emb name
%append configuration matching gui to it
%first one gets sliceoutput=true
%close it
paramfilename=[outputdirectory,'/paramfile_',handles.embryoname,suffix,num2str(i,'%04d'),'.txt'];

%paramfilename=[get(handles.outputdir,'String'),'\paramfile_',num2str(i,'%04d'),'_',suffix,'.txt'];
wrongslashes=findstr(paramfilename,'\');
paramfilename(wrongslashes)='/';

parameterfilenames{i}=paramfilename;
copyfile(get(handles.parameterfilename,'String'),paramfilename);

file=fopen(paramfilename,'a');

if((handles.lsm))
    fprintf(file,'newscope=false;\n\r');
    fprintf(file,'SIMPLETIFF=false;\n\r');
    fprintf(file,'MATLAB_STACK=false;\n\r');
    fprintf(file,'LSM_time=true;\n\r');
    
else
    if(get(handles.matlabim,'Value'))
        fprintf(file,'newscope=false;\n\r');
        fprintf(file,'SIMPLETIFF=false;\n\r');
        fprintf(file,'MATLAB_STACK=true;\n\r');
        
    else
        if(get(handles.simpletiff,'Value'))
            fprintf(file,'newscope=false;\n\r');
            fprintf(file,'MATLAB_STACK=false;\n\r');
            fprintf(file,'SIMPLETIFF=true;\n\r');

        end
    end 
end

fprintf(file,'%%Parameter overwrites generated by ROI interface:\n\r');
if(get(handles.green,'Value'))
    fprintf(file,'rednuclei=false;\n\r');
     fprintf(file,'LSM_channel=1;\n\r');
else
    fprintf(file,'rednuclei=true;\n\r');
      fprintf(file,'LSM_channel=2;\n\r');
end
fprintf(file,['start_time=',get(handles.starttime,'String'),';\n\r']);
fprintf(file,['end_time=',get(handles.endtime,'String'),';\n\r']);

if(i==1&&get(handles.makeslices,'Value'))
    fprintf(file,'outputSlice=true;\n\r');
else
    fprintf(file,'outputSlice=false;\n\r');
end

%output ROI i
fprintf(file,'ROI=true;\n\r');
points=round(handles.ROIs{1}.getPosition());
fprintf(file,['ROIxmin=',num2str(max(1,min(points(:,1)))),';\n\r']);
fprintf(file,['ROIxmax=',num2str(max(points(:,1))),';\n\r']);
fprintf(file,['ROIymin=',num2str(max(1,min(points(:,2)))),';\n\r']);
fprintf(file,['ROIymax=',num2str(max(points(:,2))),';\n\r']);
fprintf(file,['conservememory = false;\n\r']);

%build up a string for the polygonal ROI to put in parameter file
bigstring='ROIpoints=[';
for p=1:size(points,1);
    bigstring=[bigstring,num2str(points(p,1)),' ',num2str(points(p,2)),' ; '];
end
bigstring=[bigstring,']'];
fprintf(file,[bigstring,';\n\r']);%save point array in parameter file

fclose(file);
end

for i=1:length(handles.ROIs)
%for each embryo ROI
%call main detction with images 
%call tracking
%put everything in right final location
%this is done with call to central driver script
embryodir=handles.embryodir;
wrongslashes=findstr(embryodir,'\');
embryodir(wrongslashes)='/';
outputdir=get(handles.outputdir,'String');
if(iscell(outputdir))
    outputdir=outputdir{1};
end
wrongslashes=findstr(outputdir,'\');
outputdir(wrongslashes)='/';
outputdir=[outputdir,'/'];
points=round(handles.ROIs{i}.getPosition());

detect_track_driver_allmatlab_REALTIME_fish_3(parameterfilenames{i},embryodir,handles.embryoname,[suffix,num2str(i)],outputdirectory,points,false,'');%,lineageparameterfile);

end
'all embryos completed'

% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in colorchannel.
function colorchannel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in colorchannel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
%if changed update to one of default param files for red or green 
if(get(handles.green,'Value'))
set(handles.parameterfilename,'String','L:\bin\starryniteII\matlab-parameters-file-newscope_middle_weak_thresholds_journalV_besthack_nodata.txt');
else
    set(handles.parameterfilename,'String','L:\bin\starryniteII\Red_anthonytest.txt');
end


function embryosuffix_Callback(hObject, eventdata, handles)
% hObject    handle to embryosuffix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of embryosuffix as text
%        str2double(get(hObject,'String')) returns contents of embryosuffix as a double



% --- Executes during object creation, after setting all properties.
function embryosuffix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to embryosuffix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function parameterfilename_Callback(hObject, eventdata, handles)
% hObject    handle to parameterfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of parameterfilename as text
%        str2double(get(hObject,'String')) returns contents of parameterfilename as a double



% --- Executes during object creation, after setting all properties.
function parameterfilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parameterfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in parameterbutton.
function parameterbutton_Callback(hObject, eventdata, handles)
% hObject    handle to parameterbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName,FilterIndex] = uigetfile('*.*');
set(handles.parameterfilename,'String',[PathName,FileName]);



function imagefilename_Callback(hObject, eventdata, handles)
% hObject    handle to imagefilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imagefilename as text
%        str2double(get(hObject,'String')) returns contents of imagefilename as a double


 
% --- Executes during object creation, after setting all properties.
function imagefilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imagefilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in imagebutton.
function imagebutton_Callback(hObject, eventdata, handles)
% hObject    handle to imagebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName,FilterIndex] = uigetfile();
%looks for the 
ind=max(strfind(FileName,'_'));
if (isempty(ind))
    ind=regexp(FileName,'\d' ,'once');
end
ind3=(findstr(PathName,'\'));
ind3=ind3(length(ind3)-1);%take second to last occurance of dir delimater
ind2=max(findstr(FileName,'.'));

set(handles.imagefilename,'String',[PathName,FileName(1:ind-1)]);
suffix=FileName(ind2+1:length(FileName));
framenum=str2num(FileName(ind+2:ind2-1));
 handles.lsm=false;

if (strcmp(suffix,'mat'))
    load([PathName,FileName]);
    X=stack;
    clear stack;
    handles.matlabformat=true;
    im=max(X,[],3);
else
    if(strcmp(suffix,'lsm'))
        handles.lsm=true;
         handles.matlabformat=false;
         set(handles.imagefilename,'String',[PathName,FileName]);
         X=(loadCellStackLSMtime([PathName,FileName],1,1,inf));
        im=max(X,[],3);
    else
        
        if (get(handles.simpletiff,'Value'))
            X=loadSimpleStackTiff([PathName,FileName]);
            'simple tiff'
            im=max(X,[],3);
        else
            if(get(handles.green,'Value'))
                X=im2double(((loadCellStackMetamorph([PathName,FileName(1:ind-1)],framenum,2,30,[0,0,0,0],false))));
            else
                X=im2double(((loadCellStackMetamorph([PathName,FileName(1:ind-1)],framenum,1,30,[0,0,0,0],false))));
            end
            im=max(X,[],3);
            im=fliplr(im);
        end
    end
end


%im=medfilt2(im,[5,5],'symmetric');

minvalr=prctile(reshape(im,[1,numel(im)]),25);
maxvalr=prctile(reshape(im,[1,numel(im)]),99); 
im(im>maxvalr)=maxvalr;
im(im<minvalr)=minvalr;

hold on;
i=imagesc(im);
axis tight;
set(handles.axes1,'YDir','reverse');
colormap(gray);
set(i,'HitTest','off');
handles.embryodirname=PathName(ind3+1:length(PathName)-1);%end dir name
handles.embryodir=PathName;%vs full path
if(isfield(handles,'lsm')&&handles.lsm)
    handles.embryoname=FileName;
else
handles.embryoname=FileName(1:ind-1);
end
guidata(hObject,handles)


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ROIcount=handles.ROIcount+1;
handles.ROIs{handles.ROIcount}=impoly(handles.axes1);
guidata(hObject,handles)



function endtime_Callback(hObject, eventdata, handles)
% hObject    handle to endtime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endtime as text
%        str2double(get(hObject,'String')) returns contents of endtime as a double


% --- Executes during object creation, after setting all properties.
function endtime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endtime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function starttime_Callback(hObject, eventdata, handles)
% hObject    handle to starttime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of starttime as text
%        str2double(get(hObject,'String')) returns contents of starttime as a double


% --- Executes during object creation, after setting all properties.
function starttime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to starttime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in outputdirbutton.
function outputdirbutton_Callback(hObject, eventdata, handles)
% hObject    handle to outputdirbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path = uigetdir('C:\Users\shahp2\Desktop\SSD_Scratch_Space\');
set(handles.outputdir,'String',path);


function outputdir_Callback(hObject, eventdata, handles)
% hObject    handle to outputdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outputdir as text
%        str2double(get(hObject,'String')) returns contents of outputdir as a double


% --- Executes during object creation, after setting all properties.
function outputdir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outputdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in makeslices.
function makeslices_Callback(hObject, eventdata, handles)
% hObject    handle to makeslices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of makeslices



function TargetT_Callback(hObject, eventdata, handles)
% hObject    handle to TargetT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TargetT as text
%        str2double(get(hObject,'String')) returns contents of TargetT as a double


% --- Executes during object creation, after setting all properties.
function TargetT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TargetT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TargetX_Callback(hObject, eventdata, handles)
% hObject    handle to TargetX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TargetX as text
%        str2double(get(hObject,'String')) returns contents of TargetX as a double


% --- Executes during object creation, after setting all properties.
function TargetX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TargetX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TargetY_Callback(hObject, eventdata, handles)
% hObject    handle to TargetY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TargetY as text
%        str2double(get(hObject,'String')) returns contents of TargetY as a double


% --- Executes during object creation, after setting all properties.
function TargetY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TargetY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TargetZ_Callback(hObject, eventdata, handles)
% hObject    handle to TargetZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TargetZ as text
%        str2double(get(hObject,'String')) returns contents of TargetZ as a double


% --- Executes during object creation, after setting all properties.
function TargetZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TargetZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Fire.
function Fire_Callback(hObject, eventdata, handles)
% hObject    handle to Fire (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function FireDuration_Callback(hObject, eventdata, handles)
% hObject    handle to FireDuration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FireDuration as text
%        str2double(get(hObject,'String')) returns contents of FireDuration as a double


% --- Executes during object creation, after setting all properties.
function FireDuration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FireDuration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Pause.
function Pause_Callback(hObject, eventdata, handles)
% hObject    handle to Pause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Pause



function CustomAxis_Callback(hObject, eventdata, handles)
% hObject    handle to CustomAxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CustomAxis as text
%        str2double(get(hObject,'String')) returns contents of CustomAxis as a double


% --- Executes during object creation, after setting all properties.
function CustomAxis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CustomAxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AxisRebuild.
function AxisRebuild_Callback(hObject, eventdata, handles)
% hObject    handle to AxisRebuild (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
import org.rhwlab.acetree.AceTreeLauncher
import org.rhwlab.acetree.AceTree
import org.rhwlab.snight.NucleiMgr

AT = AceTree.getAceTree([]);
NM = AT.getNucleiMgr;

iIdentity = NM.getIdentity();
iAxis = iIdentity.getAxis();
iAxis = get(handles.CustomAxis,'string');
AT.buildTree(true);



function TargetCell_Callback(hObject, eventdata, handles)
% hObject    handle to TargetCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TargetCell as text
%        str2double(get(hObject,'String')) returns contents of TargetCell as a double


% --- Executes during object creation, after setting all properties.
function TargetCell_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TargetCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
