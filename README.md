Soil2D is a visualization tool for scientists to explore the effects of groundwater and soil characteristics on pollution transport.


% SOL3D M-file for Sol3D.fig

%% Programmer:  Sormeh Kashef
%      SOL3D, by itself, creates a new SOL3D or raises the existing
%      singleton*.
%
%      H = SOL3D returns the handle to a new SOL3D or the handle to
%      the existing singleton*.
%
%      SOL3D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SOL3D.M with the given input arguments.
%
%      SOL3D('Property','Value',...) creates a new SOL3D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Sol3D_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Sol3D_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Last Modified by GUIDE v2.5 29-Oct-2008 17:25:42

% Inorder to be compatible with MATLAB 7.1: In the GUIDE tool menue, GUI
% options, Resize behaviour ---> Proportional  The problem before changing
% this was the resizing of 3D image
