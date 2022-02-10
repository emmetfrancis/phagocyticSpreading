function varargout = SpreadGUI(varargin)
% SPREADGUI MATLAB code for SpreadGUI.fig
%      SPREADGUI, by itself, creates a new SPREADGUI or raises the existing
%      singleton*.
%
%      H = SPREADGUI returns the handle to a new SPREADGUI or the handle to
%      the existing singleton*.
%
%      SPREADGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPREADGUI.M with the given input arguments.
%
%      SPREADGUI('Property','Value',...) creates a new SPREADGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SpreadGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SpreadGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SpreadGUI

% Last Modified by GUIDE v2.5 02-May-2021 15:43:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SpreadGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SpreadGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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


% --- Executes just before SpreadGUI is made visible.
function SpreadGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SpreadGUI (see VARARGIN)

% Choose default command line output for SpreadGUI
handles.output = hObject;
handles.r0Adh = str2double(get(handles.r0Text,'String'));
handles.N = str2double(get(handles.NText,'String'));
handles.f0 = str2double(get(handles.f0Text,'String'));
handles.fPoly = str2double(get(handles.fPolyText,'String'));
handles.dt = str2double(get(handles.dtText,'String'));
handles.dtRef = handles.dt * .2;
handles.totalTime = str2double(get(handles.totalTimeText,'String'));
handles.numIt = round(handles.totalTime/handles.dt);
handles.mu = str2double(get(handles.muText,'String'));
handles.polyForceDecay = str2double(get(handles.polyForceDecayText,'String'));
handles.PVMult = str2double(get(handles.PVMultText,'String'));
handles.numTerms = str2double(get(handles.numTermsText,'String'));
handles.volTol = str2double(get(handles.volTolText,'String'));
handles.discreteLogic = get(handles.discreteCheckbox,'Value');
handles.PVExp = str2double(get(handles.PVExpText,'String'));
handles.timeVec = zeros(1,round(handles.totalTime/.2));
handles.PVStart = 0.01;
handles.aFactor = .001;
handles.NShift = 0;
handles.timeDecay = inf;
handles.random = false;
handles.v_rStored = {};
handles.v_zStored = {};
handles.pStored = {};
handles.meshRef = .2;
handles.meshRefineTest = false;
handles.meshVal = 0.2;
handles.epsFactor = 1;
choice = questdlg('Start with sphere or choose shape from past analysis?', 'Starting shape?',...
    'Sphere', 'Load past shape','Other test', 'Sphere');
handles.d0 = 8.5e-4;
switch choice
    case 'Sphere'
        % parameterize boundary points for half shell, in this case, want finer
        % material points for area close to the surface
        maxAdhDist = 15*handles.r0Adh*1e-7/handles.d0;
        angCrit = asin((1/2-maxAdhDist)/(1/2)) + pi/2;
        angContact = pi;
        sCrit = angCrit/2;
        sClose = angContact/2;
        sTotal = 2*pi/4;
        startSpacing = .02;
        minSpacing = handles.meshVal*handles.r0Adh*1e-7/handles.d0;
        %     L = log10(2);
        %     sCrit = (10^L * (startSpacing) + (1-10^L)*sClose) / (1-10^L);
        %     numNodes1 = round(sCrit/startSpacing);
        m = round(startSpacing/minSpacing) + 1;
        sTrans = minSpacing;
        n = 0;
        transSpan = m*(m-1)*minSpacing / 2;
        if transSpan > sClose
            transSpan = sClose;
        end
        while sTrans(end) < transSpan-minSpacing
            n = n+1;
            sTrans = [sTrans,sTrans(end)+n*minSpacing];
        end
        sTrans = sClose-sTrans;
        sTrans = sTrans(end:-1:1);
        %     numNodes2 = (sClose-sCrit)/(.2*handles.r0Adh*1e-7/handles.d0);
        %     interval1 = sCrit / (numNodes1-1);
        %     interval2 = (sClose-sCrit) / (numNodes2);
        %     sParam = [0:interval1:sCrit, sCrit+interval2:interval2:sClose];
        sStartTrans = sTrans(1);
        roundSpacing = sStartTrans/round(sStartTrans/startSpacing);
        sStartNew = 0:roundSpacing:sStartTrans-roundSpacing;
        sParam = [sStartNew, sTrans, sClose];
        %s, r, and z normalized by d0
        handles.r = .5*sin(pi*sParam/sTotal);
        %     c = 2^(1/6) - .08 / (handles.r0Adh*1e-7/(handles.d0/2));
        c = -.05;
        %s, r, and z normalized by d0
        handles.r = .5*sin(pi*sParam/sTotal);
%         c = 2^(1/6) - .08 / (handles.r0Adh*1e-7/(handles.d0/2));
        handles.z = .5*cos(pi*sParam/sTotal) + c*handles.r0Adh*1e-7/handles.d0;
        handles.z = handles.z + 0.5;% - 2^(1/6)*handles.r0Adh*1e-7/handles.d0;
        handles.z(handles.z < 0) = 0;
        handles.rStored = {};
        handles.zStored = {};
        handles.rStored{1} = handles.r;
        handles.zStored{1} = handles.z;
        handles.coeffStored = {};
        handles.noSurf = false;
        handles.squishTest = false;
        handles.squishForce = 0;
    case 'Load past shape'
        dataPath = uigetdir('','Open folder with saved data');
        rStruct = load(fullfile(dataPath,'rData.mat'));
        zStruct = load(fullfile(dataPath,'zData.mat'));
        timeStruct = load(fullfile(dataPath,'timeVec.mat'));
        handles.rStored = rStruct.rStored;
        maxIdx = length(handles.rStored);
        for i = 2:length(handles.rStored)
            if isempty(handles.rStored{i})
                maxIdx = i-1;
                break
            end
        end
        handles.rStored = handles.rStored(1:maxIdx);
        handles.zStored = zStruct.zStored(1:maxIdx);
        handles.timeVec = timeStruct.timeVec(1:maxIdx);
        try coeffStruct = load(fullfile(dataPath,'coeffVals.mat'));
            handles.coeffStored = coeffStruct.coeffStored(1:maxIdx);
        catch
            pStruct = load(fullfile(dataPath,'pData.mat'));
            v_rStruct = load(fullfile(dataPath,'v_rData.mat'));
            v_zStruct = load(fullfile(dataPath,'v_zData.mat'));
            handles.pStored = pStruct.pStored(1:maxIdx-1);
            handles.v_rStored = v_rStruct.v_rStored(1:maxIdx-1);
            handles.v_zStored = v_zStruct.v_zStored(1:maxIdx-1);
            set(handles.FEM,'Value',1)
        end
        handles.r = handles.rStored{1};
        handles.z = handles.zStored{1};
        
        % read in log file, if it's there
        try fid = fopen(fullfile(dataPath,'simSettings.txt'),'r');
            settingsCell = textscan(fid,'%s','Delimiter','\n');
            settingsCell = settingsCell{1};
            dtStr = settingsCell{2};
            startdt = strfind(dtStr,'=');
            handles.dt = str2double(dtStr(startdt+2:end-2));
            set(handles.dtText,'String',sprintf('%.2f',handles.dt))
            NStr = settingsCell{3};
            startN = strfind(NStr,'=');
            handles.N = str2double(NStr(startN+2:end));
            set(handles.NText,'String',sprintf('%.2f',handles.N))
            r0Str = settingsCell{4};
            startr0 = strfind(r0Str,'=');
            handles.r0Adh = str2double(r0Str(startr0+2:end-3));
            set(handles.r0Text,'String',sprintf('%.2f',handles.r0Adh))
            f0Str = settingsCell{5};
            startf0 = strfind(f0Str,'=');
            handles.f0 = str2double(f0Str(startf0+2:end));
            set(handles.f0Text,'String',sprintf('%.2f',handles.f0))
            fPolyStr = settingsCell{6};
            startfPoly = strfind(fPolyStr,'=');
            handles.fPoly = str2double(fPolyStr(startfPoly+2:end));
            set(handles.fPolyText,'String',sprintf('%.2f',handles.fPoly))
            muStr = settingsCell{7};
            startmu = strfind(muStr,'=');
            handles.mu = str2double(muStr(startmu+2:end));
            set(handles.muText,'String',sprintf('%.2f',handles.mu))
            polyForceDecayStr = settingsCell{8};
            startpolyForceDecay = strfind(polyForceDecayStr,'=');
            handles.polyForceDecay = str2double(polyForceDecayStr(startpolyForceDecay+2:end));
            set(handles.polyForceDecayText,'String',sprintf('%.2f',handles.polyForceDecay))
            if ~get(handles.FEM,'Value')
                PVStr = settingsCell{9};
                startPV = strfind(PVStr,'=');
                handles.PVMult = str2double(PVStr(startPV+2:end));
                set(handles.PVMultText,'String',sprintf('.%2f',handles.PVMult))
                numTermsStr = settingsCell{10};
                startnumTerms = strfind(numTermsStr,'=');
                handles.numTerms = str2double(numTermsStr(startnumTerms+2:end));
                set(handles.numTermsText,'String',sprintf('%.2f',handles.numTerms))
            end
        end
        
%         z0ContactVals = zeros(1,length(handles.zStored));
%         for i = 1:length(handles.zStored)
%             z0ContactVals(i) = min(handles.zStored{i});
%         end
%         msgbox('Fix this code')
%         z0Contact = mean(z0ContactVals);
%         handles.r0Adh = (z0Contact+.5)*handles.d0*1e7;

        set(handles.r0Text,'String',sprintf('%.2f',handles.r0Adh))
        sParam = zeros(1,length(handles.r));
        for i = 2:length(handles.r)
            sParam(i) = sParam(i-1) + sqrt((handles.r(i)-handles.r(i-1)).^2 ...
            +(handles.z(i)-handles.z(i-1)).^2);
        end
        handles.noSurf = false;
        handles.squishTest = false;
        handles.squishForce = 0;
    case 'Other test'
        testCase = questdlg('Which type of test case?','Test case','Recovery','Squish','Recovery');
        switch testCase
            case 'Recovery'
                handles.noSurf = true;
                handles.squishTest = false;
                handles.squishForce = 0;
                pipetteSize = inputdlg('Enter size of pipette in um','Pipette size');
                handles.pipetteSize = str2double(pipetteSize{1});
                rPipette = (handles.pipetteSize * 1e-4 / 2) / handles.d0;
                PMNVol = (4/3) * pi * (.5)^3; % volume stays constant, assume 8.5um spherical cell undeformed
                endsVol = (4/3) * pi * rPipette^3; % volume of hemispherical ends
                volCyl = PMNVol - endsVol; %remaining volume is cylinder
                lengthCyl = volCyl / (pi * rPipette^2);
                % parameterize boundary points for half shell
                sTotal = pi*rPipette + lengthCyl;
                numNodes = round(sTotal / .02);
                sGap = sTotal / numNodes;
                sParam = 0:sGap:sTotal;
                r = zeros(1,length(sParam));
                z = zeros(1,length(sParam));
                logic1 = sParam < pi*rPipette/2;
                r(logic1) = rPipette * sin(sParam(logic1) / rPipette);
                z(logic1) = -rPipette * (1 - cos(sParam(logic1) / rPipette));
                logic2 = (sParam >= pi*rPipette/2) & (sParam <= (pi*rPipette/2 + lengthCyl));
                r(logic2) = rPipette;
                z(logic2) = -rPipette - (sParam(logic2) - (pi*rPipette/2));
                logic3 = sParam > (pi*rPipette/2 + lengthCyl);
                r(logic3) = rPipette * (sin((sParam(logic3) - lengthCyl) / rPipette));
                z(logic3) = -rPipette - (lengthCyl - rPipette * cos((sParam(logic3)-lengthCyl) / rPipette));
                handles.r = r;
                handles.z = z - mean(z);
                handles.rStored = {};
                handles.zStored = {};
                handles.rStored{1} = handles.r;
                handles.zStored{1} = handles.z;
                handles.coeffStored = {};
            case 'Squish'
                handles.squishTest = true;
                squishForce = inputdlg('Enter squish force','Squish force');
                handles.squishForce = str2double(squishForce{1});
                handles.noSurf = false;
                offsetAng = pi/10;
                numNodes1 = round(.5*sin(offsetAng)/.02);
                interval1 = .5*sin(offsetAng) / (numNodes1-1);
                numNodes2 = round(.5*(pi-2*offsetAng)/.02);
                interval2 = (pi-2*offsetAng)/(numNodes2-1);
                handles.r = [0:interval1:.5*sin(offsetAng), .5*sin(offsetAng+interval2:interval2:pi-offsetAng),...
                    .5*sin(pi-offsetAng)-interval1:-interval1:0];
                handles.z = [.5*cos(offsetAng)*ones(1,numNodes1), .5*cos(offsetAng+interval2:interval2:pi-offsetAng),...
                    .5*cos(pi-offsetAng)*ones(1,numNodes1-1)];
                handles.z = handles.z - min(handles.z);
                handles.rStored = {};
                handles.zStored = {};
                handles.rStored{1} = handles.r;
                handles.zStored{1} = handles.z;
                handles.coeffStored = {};
        end
end
handles.curFrame = 1;
handles = updateDisplay(handles);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SpreadGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SpreadGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function simSlider_Callback(hObject, eventdata, handles)
% hObject    handle to simSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
startSlide = get(hObject,'Min');
endSlide = get(hObject,'Max');
pos = get(hObject,'Value');
num = length(handles.rStored);
handles.curFrame = round((pos/(endSlide-startSlide))*(num-1))+1;
handles.r = handles.rStored{handles.curFrame};
handles.z = handles.zStored{handles.curFrame};
handles = updateDisplay(handles);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function simSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to simSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in startSim.
function startSim_Callback(hObject, eventdata, handles)
% hObject    handle to startSim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.simPlot)
handles.NShift = 0;
% shiftGap = sqrt(1/handles.N) / 10;
% handles.NShift = 5*shiftGap;
handles.timeDecay = inf;
if ~get(handles.FEM,'Value')
    outStruct = spreadCellCalcNew(handles);
    handles.rStored = outStruct.rStored;
    handles.zStored = outStruct.zStored;
    handles.coeffStored = outStruct.coeffStored;
    handles.errorStored = outStruct.errorStored;
    handles.PVStored = outStruct.PVStored;
    handles.timeVec = outStruct.timeVec;
    handles.r = handles.rStored{handles.curFrame};
    handles.z = handles.zStored{handles.curFrame};
else
    handles.rStored = handles.rStored(1:handles.curFrame);
    handles.zStored = handles.zStored(1:handles.curFrame);
%     handles.v_rStored = handles.v_rStored(1:handles.curFrame);
%     handles.v_zStored = handles.v_zStored(1:handles.curFrame);
%     handles.pStored = handles.pStored(1:handles.curFrame);
%     handles.timeVec = handles.timeVec(1:handles.curFrame);
    outStruct = spreadCellCalcFEM(handles);
    handles.rStored = outStruct.rStored;
    handles.zStored = outStruct.zStored;
    handles.v_zStored = outStruct.v_zStored;
    handles.v_rStored = outStruct.v_rStored;
    handles.pStored = outStruct.pStored;
    handles.timeVec = outStruct.timeVec;
    handles.tensionProt = outStruct.tensionProt;
    handles.r = handles.rStored{handles.curFrame};
    handles.z = handles.zStored{handles.curFrame};
end
guidata(hObject,handles)


function r0Text_Callback(hObject, eventdata, handles)
% hObject    handle to r0Text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of r0Text as text
%        str2double(get(hObject,'String')) returns contents of r0Text as a double
if ~handles.noSurf
    handles.r0Adh = str2double(get(hObject,'String'));
    maxAdhDist = 20*handles.r0Adh*1e-7/handles.d0;
    angCrit = asin((1/2-maxAdhDist)/(1/2)) + pi/2;
    angContact = pi;
    sCrit = angCrit/2;
    sClose = angContact/2;
    sTotal = 2*pi/4;
    startSpacing = .02;
    minSpacing = handles.meshVal*handles.r0Adh*1e-7/handles.d0;
%     L = log10(2);
%     sCrit = (10^L * (startSpacing) + (1-10^L)*sClose) / (1-10^L);
%     numNodes1 = round(sCrit/startSpacing);
    m = round(startSpacing/minSpacing) + 1;
    sTrans = minSpacing;
    n = 0;
    transSpan = m*(m-1)*minSpacing / 2;
    if transSpan > sClose
        transSpan = sClose;
    end
    while sTrans(end) < transSpan-minSpacing
        n = n+1;
        sTrans = [sTrans,sTrans(end)+n*minSpacing]; 
    end
    sTrans = sClose-sTrans;
    sTrans = sTrans(end:-1:1);
%     numNodes2 = (sClose-sCrit)/(.2*handles.r0Adh*1e-7/handles.d0);
%     interval1 = sCrit / (numNodes1-1);
%     interval2 = (sClose-sCrit) / (numNodes2);
%     sParam = [0:interval1:sCrit, sCrit+interval2:interval2:sClose];
    sStartTrans = sTrans(1);
    roundSpacing = sStartTrans/round(sStartTrans/startSpacing);
    sStartNew = 0:roundSpacing:sStartTrans-roundSpacing;
    sParam = [sStartNew, sTrans, sClose];
    %s, r, and z normalized by d0
    handles.r = .5*sin(pi*sParam/sTotal);
%     c = 2^(1/6) - .08 / (handles.r0Adh*1e-7/(handles.d0/2));
    c = -.05;
    handles.z = .5*cos(pi*sParam/sTotal) + c*handles.r0Adh*1e-7/handles.d0;
    handles.z = handles.z + 0.5;% - 2^(1/6)*handles.r0Adh*1e-7/handles.d0;
    handles.rStored = {};
    handles.zStored = {};
    handles.rStored{1} = handles.r;
    handles.zStored{1} = handles.z;
    sParam = zeros(1,length(handles.r));
    for i = 2:length(handles.r)
        sParam(i) = sParam(i-1) + sqrt((handles.r(i)-handles.r(i-1)).^2 ...
            +(handles.z(i)-handles.z(i-1)).^2);
    end
    handles = updateDisplay(handles);
end
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function r0Text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r0Text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f0Text_Callback(hObject, eventdata, handles)
% hObject    handle to f0Text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f0Text as text
%        str2double(get(hObject,'String')) returns contents of f0Text as a double
handles.f0 = str2double(get(hObject,'String'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function f0Text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f0Text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fPolyText_Callback(hObject, eventdata, handles)
% hObject    handle to fPolyText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fPolyText as text
%        str2double(get(hObject,'String')) returns contents of fPolyText as a double
handles.fPoly = str2double(get(hObject,'String'));
if handles.squishTest
    handles.squishForce = str2double(get(hObject,'String'));
    handles.rStored = handles.rStored(1);
    handles.zStored = handles.zStored(1);
    handles.r = handles.rStored{1};
    handles.z = handles.zStored{1};
    handles.coeffStored = {};
    handles.curFrame = 1;
    handles = updateDisplay(handles);
end
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function fPolyText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fPolyText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dtText_Callback(hObject, eventdata, handles)
% hObject    handle to dtText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dtText as text
%        str2double(get(hObject,'String')) returns contents of dtText as a double
handles.dt = str2double(get(hObject,'String'));
handles.dtRef = 0.2 * handles.dt;
handles.numIt = round(handles.totalTime/handles.dtRef) + 1;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function dtText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dtText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function muText_Callback(hObject, eventdata, handles)
% hObject    handle to muText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of muText as text
%        str2double(get(hObject,'String')) returns contents of muText as a double
handles.mu = str2double(get(hObject,'String'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function muText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to muText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function polyForceDecayText_Callback(hObject, eventdata, handles)
% hObject    handle to polyForceDecayText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of polyForceDecayText as text
%        str2double(get(hObject,'String')) returns contents of polyForceDecayText as a double
handles.polyForceDecay = str2double(get(hObject,'String'));
if handles.noSurf
    handles.pipetteSize = handles.polyForceDecay;
    rPipette = (handles.pipetteSize * 1e-4 / 2) / handles.d0;
    PMNVol = (4/3) * pi * (.5)^3; % volume stays constant, assume 8.5um spherical cell undeformed
    endsVol = (4/3) * pi * rPipette^3; % volume of hemispherical ends
    volCyl = PMNVol - endsVol; %remaining volume is cylinder
    lengthCyl = volCyl / (pi * rPipette^2);
    % parameterize boundary points for half shell
    sTotal = pi*rPipette + lengthCyl;
    numNodes = round(sTotal / .02);
    sGap = sTotal / numNodes;
    sParam = 0:sGap:sTotal;
    r = zeros(1,length(sParam));
    z = zeros(1,length(sParam));
    logic1 = sParam < pi*rPipette/2;
    r(logic1) = rPipette * sin(sParam(logic1) / rPipette);
    z(logic1) = -rPipette * (1 - cos(sParam(logic1) / rPipette));
    logic2 = (sParam >= pi*rPipette/2) & (sParam <= (pi*rPipette/2 + lengthCyl));
    r(logic2) = rPipette;
    z(logic2) = -rPipette - (sParam(logic2) - (pi*rPipette/2));
    logic3 = sParam > (pi*rPipette/2 + lengthCyl);
    r(logic3) = rPipette * (sin((sParam(logic3) - lengthCyl) / rPipette));
    z(logic3) = -rPipette - (lengthCyl - rPipette * cos((sParam(logic3)-lengthCyl) / rPipette));
    handles.r = r;
    handles.z = z - mean(z);
    handles.rStored = {};
    handles.zStored = {};
    handles.rStored{1} = handles.r;
    handles.zStored{1} = handles.z;
    handles.coeffStored = {};
    handles.curFrame = 1;
    handles = updateDisplay(handles);
end
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function polyForceDecayText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to polyForceDecayText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NText_Callback(hObject, eventdata, handles)
% hObject    handle to NText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NText as text
%        str2double(get(hObject,'String')) returns contents of NText as a double
handles.N = str2double(get(hObject,'String'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function NText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function totalTimeText_Callback(hObject, eventdata, handles)
% hObject    handle to totalTimeText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of totalTimeText as text
%        str2double(get(hObject,'String')) returns contents of totalTimeText as a double
handles.totalTime = str2num(get(handles.totalTimeText,'String'));
handles.numIt = round(handles.totalTime/handles.dtRef) + 1;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function totalTimeText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to totalTimeText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in saveData.
function saveData_Callback(hObject, eventdata, handles)
% hObject    handle to saveData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
savePath = uigetdir('','Choose folder to save r and z data');
if ~get(handles.FEM,'Value')
    rStored = handles.rStored;
    zStored = handles.zStored;
    PVStored = handles.PVStored;
    timeVec = handles.timeVec;
    coeffStored = handles.coeffStored;
    errorStored = handles.errorStored;
    save(fullfile(savePath,'rData.mat'),'rStored')
    save(fullfile(savePath,'zData.mat'),'zStored')
    save(fullfile(savePath,'coeffVals.mat'),'coeffStored')
    save(fullfile(savePath,'errorStored.mat'),'errorStored')
    if ~handles.noSurf && ~handles.squishTest
        save(fullfile(savePath,'PVStored.mat'),'PVStored')
        save(fullfile(savePath,'timeVec.mat'),'timeVec')
    end
    saveLog(handles,savePath)
else
    rStored = handles.rStored;
    zStored = handles.zStored;
    v_rStored = handles.v_rStored;
    v_zStored = handles.v_zStored;
    pStored = handles.pStored;
    timeVec = handles.timeVec;
    tensionProt = handles.tensionProt;
    save(fullfile(savePath,'rData.mat'),'rStored')
    save(fullfile(savePath,'zData.mat'),'zStored')
    save(fullfile(savePath,'v_rData.mat'),'v_rStored')
    save(fullfile(savePath,'v_zData.mat'),'v_zStored')
    save(fullfile(savePath,'pData.mat'),'pStored')
    save(fullfile(savePath,'timeVec.mat'),'timeVec')
    save(fullfile(savePath,'tensionProt.mat'),'tensionProt')
    saveFEMLog(handles,savePath)
end
guidata(hObject,handles)


% --- Executes on button press in batchRun.
function batchRun_Callback(hObject, eventdata, handles)
% hObject    handle to batchRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.noSurf
    testValsCell = inputdlg('Enter pipette sizes to test, separated by commas','Enter values');
    testVals = textscan(testValsCell{1},'%.6f,');
    testVals = double(testVals{:});
    savePath = uigetdir('','Choose location to save folders with batch data');
    for i = 1:length(testVals)
        handles.pipetteSize = testVals(i);
        rPipette = (handles.pipetteSize * 1e-4 / 2) / handles.d0;
        PMNVol = (4/3) * pi * (.5)^3; % volume stays constant, assume 8.5um spherical cell undeformed
        endsVol = (4/3) * pi * rPipette^3; % volume of hemispherical ends
        volCyl = PMNVol - endsVol; %remaining volume is cylinder
        lengthCyl = volCyl / (pi * rPipette^2);
        % parameterize boundary points for half shell
        sTotal = pi*rPipette + lengthCyl;
        numNodes = round(sTotal / .02);
        sGap = sTotal / numNodes;
        sParam = 0:sGap:sTotal;
        r = zeros(1,length(sParam));
        z = zeros(1,length(sParam));
        logic1 = sParam < pi*rPipette/2;
        r(logic1) = rPipette * sin(sParam(logic1) / rPipette);
        z(logic1) = -rPipette * (1 - cos(sParam(logic1) / rPipette));
        logic2 = (sParam >= pi*rPipette/2) & (sParam <= (pi*rPipette/2 + lengthCyl));
        r(logic2) = rPipette;
        z(logic2) = -rPipette - (sParam(logic2) - (pi*rPipette/2));
        logic3 = sParam > (pi*rPipette/2 + lengthCyl);
        r(logic3) = rPipette * (sin((sParam(logic3) - lengthCyl) / rPipette));
        z(logic3) = -rPipette - (lengthCyl - rPipette * cos((sParam(logic3)-lengthCyl) / rPipette));
        handles.r = r;
        handles.z = z - mean(z);
        handles.rStored = {};
        handles.zStored = {};
        handles.rStored{1} = handles.r;
        handles.zStored{1} = handles.z;
        handles.coeffStored = {};
        outStruct = spreadCellCalcNew(handles);
        rStored = outStruct.rStored;
        zStored = outStruct.zStored;
        coeffStored = outStruct.coeffStored;
        errorStored = outStruct.errorStored;
        timeVec = outStruct.timeVec;
        curFolder = sprintf('%.2fumPipette',handles.pipetteSize);
        mkdir(savePath,curFolder);
        save(fullfile(savePath,curFolder,'rData.mat'),'rStored')
        save(fullfile(savePath,curFolder,'zData.mat'),'zStored')
        save(fullfile(savePath,curFolder,'coeffVals.mat'),'coeffStored')
        save(fullfile(savePath,curFolder,'errorStored.mat'),'errorStored')
        save(fullfile(savePath,curFolder,'timeVec.mat'),'timeVec')
        saveLog(handles,fullfile(savePath,curFolder))
    end
elseif handles.squishTest
    testValsCell = inputdlg('Enter squishing forces to test, separated by commas','Enter values');
    testVals = textscan(testValsCell{1},'%.6f,');
    testVals = double(testVals{:});
    savePath = uigetdir('','Choose location to save folders with batch data');
    for i = 1:length(testVals)
        handles.squishForce = testVals(i);
        handles.rStored = handles.rStored(1);
        handles.zStored = handles.zStored(1);
        handles.r = handles.rStored{1};
        handles.z = handles.zStored{1};
        handles.coeffStored = {};
        outStruct = spreadCellCalcNew(handles);
        rStored = outStruct.rStored;
        zStored = outStruct.zStored;
        coeffStored = outStruct.coeffStored;
        errorStored = outStruct.errorStored;
        curFolder = sprintf('%.2fsquishForce',handles.squishForce);
        mkdir(savePath,curFolder);
        save(fullfile(savePath,curFolder,'rData.mat'),'rStored')
        save(fullfile(savePath,curFolder,'zData.mat'),'zStored')
        save(fullfile(savePath,curFolder,'coeffVals.mat'),'coeffStored')
        save(fullfile(savePath,curFolder,'errorStored.mat'),'errorStored')
        saveLog(handles,fullfile(savePath,curFolder))
    end
else
    handles.NShift = 0;
    handles.timeDecay = inf;
    singDoub = questdlg('Single or double batch?','Double batch?','Single','Double','Double with shift','Single');
    NShiftRepeats = 0;
    NShiftLogic = false;
    if strcmp(singDoub,'Double')
        numLoops = 2;
    elseif strcmp(singDoub,'Double with shift')
        numLoops = 2;
        shiftText = inputdlg('Enter number of shifts','Enter number of shifts');
        NShiftRepeats = str2double(shiftText{1});
        NShiftLogic = true;
    else
        numLoops = 1;
    end
    storeTestVals = cell(1,2);
    paramChoice = cell(1,2);
    for loop = 1:numLoops
        paramChoice{loop} = questdlg('Which parameter do you want to vary?','Which parameter?',...
            'N','Numerical','fPoly','N');
        if strcmp(paramChoice{loop},'NShift')
            shiftText = inputdlg('Enter number of shifts','Enter number of shifts');
            shiftVal = str2double(shiftText{1});
            testValsCell = inputdlg('Now enter N vals to test, separated by commas','Enter N vals');
            testVals = textscan(testValsCell{1},'%.6f,');
            storeTestVals{2} = double(testVals{:});
            minN = min(storeTestVals{2});
            paramChoice{2} = 'N';
            paramChoice{1} = 'NShift';
            %         shiftVal = shiftVal * 1e-3;
            %         storeTestVals{1} = 0:shiftVal:sqrt(1/minN);
            numLoops = 2;
            break
        elseif strcmp(paramChoice{loop},'Numerical')
%             dtValsCell = inputdlg('Now enter dt vals to test, separated by commas','Enter dt vals');
%             dtVals = textscan(dtValsCell{1},'%.6f,');
%             dtVals = double(dtVals{:});
            dtVals = [];
            meshValsCell = inputdlg('Now enter eps vals to test, separated by commas','Enter eps vals');
            meshVals = textscan(meshValsCell{1},'%.6f,');
            meshVals = double(meshVals{:});
            paramChoice{1} = 'Numerical';
            storeTestVals{1} = [dtVals, meshVals];
            paramChoice{2} = 'Numerical BZ vs PZ';
            storeTestVals{2} = {'BZ','PZ'};
            numLoops = 2;
        else
            testValsCell = inputdlg(sprintf('Enter values to test for %s, separated by commas',paramChoice{loop}),'Enter values');
            testVals = textscan(testValsCell{1},'%.6f,');
            storeTestVals{loop} = double(testVals{:});
        end
    end
    if numLoops == 1
        storeTestVals{2} = handles.N;
        paramChoice{2} = 'N';
    end
    if NShiftRepeats == 0
        repeats = inputdlg('How many repeats?:');
        repeats = str2double(repeats{1});
        if repeats > 1
            handles.random = true;
        else
            handles.random = false;
        end
    else
        handles.random = false;
        repeats = NShiftRepeats;
    end
    savePath = uigetdir('','Choose location to save folders with batch data');
    handles.rStart = handles.r;
    handles.zStart = handles.z;
    handles.meshRefineTest = false;
    for k = 1:repeats
        for i = 1:length(storeTestVals{2}) % if only one batch, this doesn't affect inner loop
            switch paramChoice{2}
                case 'N'
                    handles.N = storeTestVals{2}(i);
                    handles.dtRef = str2double(get(handles.dtText,'String'));
                    handles.numIt = round(handles.totalTime/handles.dtRef) + 1;
                    if strcmp(paramChoice{1},'NShift')
                        shiftGap = sqrt(1/handles.N) / shiftVal;
                        storeTestVals{1} = shiftGap:shiftGap:sqrt(1/handles.N);
                    end
                    if NShiftLogic
                        shiftGap = sqrt(1/handles.N) / NShiftRepeats;
                        handles.NShift = (k-1) * shiftGap;
                    end
                case 'r0Adh'
                    handles.r0Adh = storeTestVals{2}(i);
                case 'numTerms'
                    handles.numTerms = round(storeTestVals{2}(i));
                case 'fPoly'
                    handles.fPoly = storeTestVals{2}(i);
                case 'f0'
                    handles.f0 = storeTestVals{2}(i);
                case 'timeDecay'
                    handles.timeDecay = storeTestVals{2}(i);
                case 'Length decay'
                    handles.polyForceDecay = storeTestVals{2}(i);
                    handles.fPoly = 2400 / handles.polyForceDecay;
                case 'dt'
                    handles.dt = storeTestVals{2}(i);
                case 'Mesh Refinement'
                    handles.meshRefineTest = true;
                    handles.meshVal = storeTestVals{2}(i);
                case 'Eps factor'
                    handles.epsFactor = storeTestVals{2}(i);
                case 'Numerical BZ vs PZ'
                    switch storeTestVals{2}{i}
                        case 'BZ'
                            handles.N = 1000;
                            handles.fPoly = 0;
                            handles.mu = 2000;
                            handles.totalTime = 200;
                        case 'PZ'
                            handles.N = 0;
                            handles.fPoly = 3000;
                            handles.mu = 16600;
                            handles.totalTime = 150;
                    end
                    handles.numIt = round(handles.totalTime/handles.dtRef) + 1;
            end
            for j = 1:length(storeTestVals{1})
                handles.rStored = {};
                handles.zStored = {};
                handles.rStored{1} = handles.rStart;
                handles.zStored{1} = handles.zStart;
                handles.r = handles.rStart;
                handles.z = handles.zStart;
                switch paramChoice{1}
                    case 'N'
                        handles.N = storeTestVals{1}(j);
                        handles.dtRef = str2double(get(handles.dtText,'String'));
                        handles.numIt = round(handles.totalTime/handles.dtRef) + 1;
                        if NShiftLogic
                            shiftGap = sqrt(1/handles.N) / NShiftRepeats;
                            handles.NShift = (k-1) * shiftGap;
                        end
                    case 'r0Adh'
                        handles.r0Adh = storeTestVals{1}(j);
                    case 'numTerms'
                        handles.numTerms = round(storeTestVals{1}(j));
                    case 'NShift'
                        maxShift = sqrt(1/handles.N);
                        if storeTestVals{1}(j) < maxShift
                            handles.NShift = storeTestVals{1}(j);
                        else
                            continue;
                        end
                    case 'fPoly'
                        handles.fPoly = storeTestVals{1}(j);
                    case 'f0'
                        handles.f0 = storeTestVals{1}(j);
                    case 'timeDecay'
                        handles.timeDecay = storeTestVals{1}(j);
                    case 'Length decay'
                        handles.polyForceDecay = storeTestVals{1}(j);
                        handles.fPoly = 2400 / handles.polyForceDecay;
                    case 'dt'
                        handles.dt = storeTestVals{1}(j);
                    case 'Mesh Refinement'
                        handles.meshRefineTest = true;
                        handles.meshVal = storeTestVals{1}(j);
                    case 'Eps factor'
                        handles.epsFactor = storeTestVals{1}(j);
                    case 'Numerical'
                        if j <= length(dtVals)
                            handles.dt = storeTestVals{1}(j);
                            handles.meshRefineTest = false;
                            handles.meshVal = 0.05;
                        else
%                             handles.dt = 0.1;
%                             handles.meshRefineTest = true;
                            handles.epsFactor = storeTestVals{1}(j);
                        end
                end
                if ~get(handles.FEM,'Value')
                    outStruct = spreadCellCalcNew(handles);
                    rStored = outStruct.rStored;
                    zStored = outStruct.zStored;
                    coeffStored = outStruct.coeffStored;
                    errorStored = outStruct.errorStored;
                    PVStored = outStruct.PVStored;
                    timeVec = outStruct.timeVec;
                else
                    outStruct = spreadCellCalcFEM(handles);
                    rStored = outStruct.rStored;
                    zStored = outStruct.zStored;
                    v_zStored = outStruct.v_zStored;
                    v_rStored = outStruct.v_rStored;
                    pStored = outStruct.pStored;
                    timeVec = outStruct.timeVec;
                    tensionProt = outStruct.tensionProt;
                end
                if numLoops == 1
                    curFolder = sprintf('%s%.5f',paramChoice{1},storeTestVals{1}(j));
                else
                    if strcmp(paramChoice{1},'Numerical')
                        if j <= length(dtVals)
                            curFolder = sprintf('dt%.5f_Numerical%s',storeTestVals{1}(j),...
                                storeTestVals{2}{i});
                        else
                            curFolder = sprintf('MeshRefine%.5f_Numerical%s',storeTestVals{1}(j),...
                                storeTestVals{2}{i});
                        end
                    else
                        curFolder = sprintf('%s%.5f_%s%.5f',paramChoice{1},storeTestVals{1}(j),...
                            paramChoice{2},storeTestVals{2}(i));
                    end
                end
                if repeats > 1 && ~NShiftLogic
                    curFolder = sprintf('%s_repeat%d',curFolder,k);
                    NLoc = outStruct.NLoc;
                    mkdir(savePath,curFolder);
                    save(fullfile(savePath,curFolder,'NLoc.mat'),'NLoc')
                elseif NShiftLogic
                    curFolder = sprintf('%s_shift%dof%d',curFolder,k,NShiftRepeats);
                    mkdir(savePath,curFolder);
                else
                    mkdir(savePath,curFolder);
                end
                if ~get(handles.FEM,'Value')
                    save(fullfile(savePath,curFolder,'rData.mat'),'rStored')
                    save(fullfile(savePath,curFolder,'zData.mat'),'zStored')
                    save(fullfile(savePath,curFolder,'coeffVals.mat'),'coeffStored')
                    save(fullfile(savePath,curFolder,'PVStored.mat'),'PVStored')
                    save(fullfile(savePath,curFolder,'errorStored.mat'),'errorStored')
                    save(fullfile(savePath,curFolder,'timeVec.mat'),'timeVec')
                    saveLog(handles,fullfile(savePath,curFolder))
                else
                    save(fullfile(savePath,curFolder,'rData.mat'),'rStored')
                    save(fullfile(savePath,curFolder,'zData.mat'),'zStored')
                    save(fullfile(savePath,curFolder,'v_rData.mat'),'v_rStored')
                    save(fullfile(savePath,curFolder,'v_zData.mat'),'v_zStored')
                    save(fullfile(savePath,curFolder,'pData.mat'),'pStored')
                    save(fullfile(savePath,curFolder,'timeVec.mat'),'timeVec')
                    save(fullfile(savePath,curFolder,'tensionProt.mat'),'tensionProt')
                    saveFEMLog(handles,fullfile(savePath,curFolder))
                end
            end
        end
    end
end
handles.rStored = rStored;
handles.zStored = zStored;
handles.r = handles.rStored{handles.curFrame};
handles.z = handles.zStored{handles.curFrame};
handles.meshRefineTest = false;
guidata(hObject,handles)

% --- Executes on button press in numBatch.
function numBatch_Callback(hObject, eventdata, handles) % only enabled for cortical shell model currently
% hObject    handle to numBatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% choice = questdlg('Which parameter do you want to vary?','Which parameter?',...
%     'dt','Vol Spring','Number of Terms','dt');
singDoub = questdlg('Single or double batch?','Double batch?','Single','Double','Single');
if strcmp(singDoub,'Double')
    numLoops = 2;
else
    numLoops = 1;
end
storeTestVals = cell(1,2);
paramChoice = cell(1,2);
for loop = 1:numLoops
    paramChoice{loop} = questdlg('Which parameter do you want to vary? (numerical first, then N)','Which parameter?',...
        'N','numTerms','PVMult','N');
    testValsCell = inputdlg(sprintf('Enter values to test for %s, separated by commas',paramChoice{loop}),'Enter values');
    testVals = textscan(testValsCell{1},'%.6f,');
    storeTestVals{loop} = double(testVals{:});
end
if numLoops == 1
    storeTestVals{2} = handles.N;
    paramChoice{2} = 'N';
end
% repeats = inputdlg('How many repeats?:');
% repeats = str2double(repeats{1});
savePath = uigetdir('','Choose location to save folders with batch data');
handles.rStart = handles.r;
handles.zStart = handles.z;
% for k = 1:repeats
for i = 1:length(storeTestVals{2}) % if only one batch, this doesn't affect inner loop
    switch paramChoice{2}
        case 'N'
            handles.N = storeTestVals{2}(i);
            dtRef = str2double(get(handles.dtText,'String'));
            if handles.fPoly > 0
                handles.dt = dtRef;
            else
                if handles.N < 200
                    handles.dtRef = 2;
                    handles.PVMult = 1;
                elseif handles.N < 500
                    handles.dtRef = 1;
                    handles.PVMult = 1;
                elseif handles.N < 1500
                    handles.dtRef = 0.5;
                    handles.PVMult = 1;
                elseif handles.N < 5000
                    handles.dtRef = 0.1;
                    handles.PVMult = 1;
                else
                    handles.dtRef = .02;
                    handles.PVMult = 1;
                end
%                 handles.dt = dtRef*(2*log((max(storeTestVals{1})/handles.N))+1);
            end
            handles.numIt = round(handles.totalTime/handles.dtRef) + 1;
        case 'PVExp'
            handles.PVExp = storeTestVals{2}(i);
        case 'PVMult'
            handles.PVMult = storeTestVals{2}(i);
        case 'PVStart'
            handles.PVStart = storeTestVals{2}(i);
        case 'dt'
            handles.dt = storeTestVals{2}(i);
            handles.dtRef = handles.dtRef * handles.dt;
            handles.numIt = round(handles.totalTime/handles.dtRef) + 1;
        case 'numTerms'
            handles.numTerms = round(storeTestVals{2}(i));
        case 'aFactor'
            handles.aFactor = storeTestVals{2}(i);
    end
    for j = 1:length(storeTestVals{1})
        handles.rStored = {};
        handles.zStored = {};
        handles.coeffStored = {};
        handles.rStored{1} = handles.rStart;
        handles.zStored{1} = handles.zStart;
        handles.r = handles.rStart;
        handles.z = handles.zStart;
        switch paramChoice{1}
            case 'N'
                handles.N = storeTestVals{1}(j);
                dtRef = str2double(get(handles.dtText,'String'));
                if handles.fPoly > 0
                    handles.dt = dtRef;
                else
%                     handles.dt = dtRef*(2*log((max(storeTestVals{1})/handles.N))+1);
                    if handles.N < 200
                        handles.dtRef = 2;
                        handles.PVMult = 1;
                    elseif handles.N < 500
                        handles.dtRef = 1;
                        handles.PVMult = 1;
                    elseif handles.N < 1500
                        handles.dtRef = 0.5;
                        handles.PVMult = 1;
                    elseif handles.N < 5000
                        handles.dtRef = 0.1;
                        handles.PVMult = 1;
                    else
                        handles.dtRef = 0.02;
                        handles.PVMult = 1;
                    end
                end
                handles.numIt = round(handles.totalTime/handles.dtRef) + 1;
            case 'PVExp'
                handles.PVExp = storeTestVals{1}(j);
            case 'PVMult'
                handles.PVMult = storeTestVals{1}(j);
            case 'PVStart'
                handles.PVStart = storeTestVals{1}(j);
            case 'dt'
                handles.dt = storeTestVals{1}(j);
                handles.dtRef = handles.dtRef * handles.dt;
                handles.numIt = round(handles.totalTime/handles.dtRef) + 1;
            case 'numTerms'
                handles.numTerms = round(storeTestVals{1}(j));
            case 'aFactor'
                handles.aFactor = storeTestVals{1}(j);
        end
        if ~get(handles.FEM,'Value')
            outStruct = spreadCellCalcNew(handles);
            rStored = outStruct.rStored;
            zStored = outStruct.zStored;
            coeffStored = outStruct.coeffStored;
            errorStored = outStruct.errorStored;
            PVStored = outStruct.PVStored;
            timeVec = outStruct.timeVec;
        else
            outStruct = spreadCellCalcFEM(handles);
            rStored = outStruct.rStored;
            zStored = outStruct.zStored;
            v_zStored = outStruct.v_zStored;
            v_rStored = outStruct.v_rStored;
            pStored = outStruct.pStored;
            timeVec = outStruct.timeVec;
        end
        if numLoops == 1
            curFolder = sprintf('%s%.3f',paramChoice{1},storeTestVals{1}(j));
            if strcmp(paramChoice{1}, 'dt')
                curFolder = sprintf('%s%.3fMult',paramChoice{1},storeTestVals{1}(j));
            end
        else
            curFolder = sprintf('%s%.3f_%s%.3f',paramChoice{1},storeTestVals{1}(j),...
                paramChoice{2},storeTestVals{2}(i));
            if strcmp(paramChoice{1}, 'dt')
                curFolder = sprintf('%s%.3fMult_%s%.3f',paramChoice{1},storeTestVals{1}(j),...
                paramChoice{2},storeTestVals{2}(i));
            end
        end
        mkdir(savePath,curFolder);
        if ~get(handles.FEM,'Value')
            save(fullfile(savePath,curFolder,'rData.mat'),'rStored')
            save(fullfile(savePath,curFolder,'zData.mat'),'zStored')
            save(fullfile(savePath,curFolder,'coeffVals.mat'),'coeffStored')
            save(fullfile(savePath,curFolder,'PVStored.mat'),'PVStored')
            save(fullfile(savePath,curFolder,'errorStored.mat'),'errorStored')
            save(fullfile(savePath,curFolder,'timeVec.mat'),'timeVec')
            saveLog(handles,fullfile(savePath,curFolder))
        else
            save(fullfile(savePath,curFolder,'rData.mat'),'rStored')
            save(fullfile(savePath,curFolder,'zData.mat'),'zStored')
            save(fullfile(savePath,curFolder,'v_rData.mat'),'v_rStored')
            save(fullfile(savePath,curFolder,'v_zData.mat'),'v_zStored')
            save(fullfile(savePath,curFolder,'pData.mat'),'pStored')
            save(fullfile(savePath,curFolder,'timeVec.mat'),'timeVec')
            saveFEMLog(handles,fullfile(savePath,curFolder))
        end
    end
end
% end
handles.rStored = rStored;
handles.zStored = zStored;
handles.r = handles.rStored{handles.curFrame};
handles.z = handles.zStored{handles.curFrame};
guidata(hObject,handles)


function PVMultText_Callback(hObject, eventdata, handles)
% hObject    handle to PVMultText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PVMultText as text
%        str2double(get(hObject,'String')) returns contents of PVMultText as a double
handles.PVMult = str2double(get(handles.PVMultText,'String'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function PVMultText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PVMultText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function numTermsText_Callback(hObject, eventdata, handles)
% hObject    handle to numTermsText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numTermsText as text
%        str2double(get(hObject,'String')) returns contents of numTermsText as a double
handles.numTerms = str2num(get(hObject,'String'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function numTermsText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numTermsText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function handles = updateDisplay(handles)
    axes(handles.simPlot)
    cla
    hold on
    title(sprintf('%.2fs',handles.timeVec(handles.curFrame)))
    r = handles.r;
    z = handles.z;
%     plot([-1.5 1.5], [-.5 -.5], 'k')
    if handles.noSurf
        plot(z, r, 'r*')
        plot(z, -r, 'r*')
        xlim([-1.5 1.5])
        ylim([-.6 .6])
    elseif handles.squishTest
        plot(r, z, 'r*')
        plot(-r, z, 'r*')
        xLower = [-1.5, -1.5, 1.5, 1.5];
        yLower = [-.1, 0, 0, -.1];
        patch('XData',xLower,'YData',yLower,'FaceColor',[.5 .5 .5],'EdgeColor','k')
        xUpper = [-1.5, -1.5, 1.5, 1.5];
        yUpper = [max(z), max(z)+.1, max(z)+.1, max(z)];
        patch('XData',xUpper,'YData',yUpper,'FaceColor',[.5 .5 .5],'EdgeColor','k')
        xlim([-1.5 1.5])
        ylim([-.1 1.5])
    else
        if ~get(handles.FEM,'Value')
            plot(r, z, 'r*')
            plot(-r, z, 'r*')
        else
            numCols = sum(z <= eps);
            numRows = length(r)/numCols;
            xMesh = reshape(r,[numCols,numRows])';
            yMesh = reshape(z,[numCols,numRows])';
            surf(xMesh,yMesh,zeros(size(xMesh)),'FaceColor','none','EdgeColor','k')
            surf(-xMesh,yMesh,zeros(size(xMesh)),'FaceColor','none','EdgeColor','k')
            view([0 90])
        end
        x = [-1.5, -1.5, 1.5, 1.5];
        y = [-.1, 0, 0, -.1];
        patch('XData',x,'YData',y','FaceColor','k')
        xlim([-1.5 1.5])
        ylim([-.1 1.5])
    end
    daspect([1 1 1])
    drawnow

function saveFEMLog(handles,savePath)
    outString{1} = sprintf('Number of iterations: %d', handles.numIt);
    outString{2} = sprintf('dt = %.4f s',handles.dt);
    outString{3} = sprintf('N = %.2f',handles.N);
    outString{4} = sprintf('r0 = %.2f nm',handles.r0Adh);
    outString{5} = sprintf('f0 = %.4f',handles.f0);
    outString{6} = sprintf('fPoly = %.2f',handles.fPoly);
    outString{7} = sprintf('mu (viscosity) = %.2f',handles.mu);
    outString{8} = sprintf('Poly force decay (um) = %.2f', handles.polyForceDecay);
    outString{9} = sprintf('FEM Solution');
    fid = fopen(fullfile(savePath,'simSettings.txt'),'w');
    fprintf(fid,'%s\r\n',outString{:});
    fclose(fid);


function saveLog(handles,savePath)
    if handles.noSurf
        outString{1} = sprintf('Number of iterations: %d', handles.numIt);
        outString{2} = sprintf('dt = %.4f s',handles.dt);
        outString{3} = sprintf('Pipette size = %.4f um', handles.pipetteSize);
        outString{4} = sprintf('Number of Terms = %d',handles.numTerms);
    elseif handles.squishTest
        outString{1} = sprintf('Number of iterations: %d', handles.numIt);
        outString{2} = sprintf('dt = %.4f s',handles.dt);
        outString{3} = sprintf('Relative squish force = %.4f um', handles.squishForce);
        outString{4} = sprintf('Number of Terms = %d',handles.numTerms);
        outString{5} = sprintf('PVMult = %.4f',handles.PVMult);
    else
        outString{1} = sprintf('Number of iterations: %d', handles.numIt);
        outString{2} = sprintf('dt = %.4f s',handles.dt);
        outString{3} = sprintf('N = %.2f',handles.N);
        outString{4} = sprintf('r0 = %.2f nm',handles.r0Adh);
        outString{5} = sprintf('f0 = %.4f',handles.f0);
        outString{6} = sprintf('fPoly = %.2f',handles.fPoly);
        outString{7} = sprintf('mu (viscosity) = %.2f',handles.mu);
        outString{8} = sprintf('Poly force decay (um) = %.2f', handles.polyForceDecay);
        outString{9} = sprintf('PVMult = %.4f',handles.PVMult);
        outString{10} = sprintf('Number of Terms = %d',handles.numTerms);
    end
    fid = fopen(fullfile(savePath,'simSettings.txt'),'w');
    fprintf(fid,'%s\r\n',outString{:});
    fclose(fid);

function volTolText_Callback(hObject, eventdata, handles)

handles.volTol = str2double(get(hObject,'String'));
guidata(hObject,handles);

function volTolText_CreateFcn(hObject, eventdata, handles)

function discreteCheckbox_Callback(hObject, eventdata, handles)

handles.discreteLogic = get(hObject,'Value');
guidata(hObject,handles)

function PVExpText_Callback(hObject, eventdata, handles)
% hObject    handle to PVExpText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PVExpText as text
%        str2double(get(hObject,'String')) returns contents of PVExpText as a double
handles.PVExp = str2double(get(hObject,'String'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function PVExpText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PVExpText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in protTest.
function protTest_Callback(hObject, eventdata, handles)
% hObject    handle to protTest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [rFile,rPath] = uigetfile('*.mat','Load reference r data');
    [zFile,zPath] = uigetfile('*.mat','Load reference z data');
    [coeffFile,coeffPath] = uigetfile('*.mat','Load reference coeff data');
    rStruct = load(fullfile(rPath,rFile));
    zStruct = load(fullfile(zPath,zFile));
    coeffStruct = load(fullfile(coeffPath,coeffFile));
    handles.rStored = rStruct.rStored;
    handles.zStored = zStruct.zStored;
    handles.coeffStored = coeffStruct.coeffStored;
    startPosText = inputdlg('Enter the starting position (in microns)');
    startPos = str2double(startPosText{1})/8.5;
    closestDist = 0.1;
    for i = 1:length(handles.rStored)
        rCur = handles.rStored{i};
        zCur = handles.zStored{i};
        if isempty(rCur)
            break
        end
        curRad = rCur(find(zCur <= 0, 1, 'first'));
        if abs(curRad-startPos) < closestDist
            closestDist = abs(curRad-startPos);
            closeIdx = i;
        end
    end
    firstContact = find(handles.zStored{closeIdx} <= 0, 1, 'first');
    handles.rStored{closeIdx}(firstContact) = startPos;
    handles.rStored = handles.rStored(closeIdx);
    handles.zStored = handles.zStored(closeIdx);
    handles.r = handles.rStored{1};
    handles.z = handles.zStored{1};
    handles.curFrame = 1;
    handles = updateDisplay(handles);
    NValsCell = inputdlg('Now enter the N values to be tested (in #/um^2), separated by commas');
    NVals = textscan(NValsCell{1},'%.6f,');
    NVals = double(NVals{:});
    handles.totalTime = 60;
    handles.numIt = round(handles.totalTime/handles.dt);
    savePath = uigetdir('','Choose parent folder for saving');
    for i = 1:length(NVals)
        handles.N = NVals(i);
        handles.NShift = mod(startPos*8.5,sqrt(1/handles.N));
        outStruct = spreadCellCalcNew(handles);
        rStored = outStruct.rStored;
        zStored = outStruct.zStored;
        timeVec = outStruct.timeVec;
        coeffStored = outStruct.coeffStored;
        PVStored = outStruct.PVStored;
        errorStored = outStruct.errorStored;
        curFolder = sprintf('start%.2fum_N%d',startPos*8.5,handles.N);
        mkdir(savePath,curFolder);
        save(fullfile(savePath,curFolder,'rData.mat'),'rStored')
        save(fullfile(savePath,curFolder,'zData.mat'),'zStored')
        save(fullfile(savePath,curFolder,'coeffVals.mat'),'coeffStored')
        save(fullfile(savePath,curFolder,'PVStored.mat'),'PVStored')
        save(fullfile(savePath,curFolder,'errorStored.mat'),'errorStored')
        save(fullfile(savePath,curFolder,'timeVec.mat'),'timeVec')
        saveLog(handles,fullfile(savePath,curFolder))
    end
    handles.rStored = rStored;
    handles.zStored = zStored;
    handles.r = handles.rStored{handles.curFrame};
    handles.z = handles.zStored{handles.curFrame};
    guidata(hObject,handles)


% --- Executes on button press in FEM.
function FEM_Callback(hObject, eventdata, handles)
% hObject    handle to FEM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FEM



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
