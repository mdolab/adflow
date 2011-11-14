%
%
% general, preliminary polara plot
%

clear all;close all;
d2r=pi/180;

%------------------------------
%
% prepare polara file
%
[status,res]=system('polaraPlotPrep.py');
if status ==0
    helpdlg(res);
else
    helpdlg('UUPS !!! polaraPlotPrep.py not found in this system. Let us hope that the polara file already exists');
end
%
% read input
%


filename=uigetfile('*PolaraCoeff.dat','Polara output file ');
if size(filename,2) < 2
    clear all;
    error('Polara file name was not specified. Think again ');
end

trailer=filename(end-14:end);
if trailer == 'PolaraCoeff.dat'
    Confcase=filename(1:end-15);
    if Confcase(end)=='_',Confcase=Confcase(1:end-1);end;
else
    Confcase=filename;
end
ii=find(Confcase=='_');
if isempty(ii)==0,Confcase(ii)='-';end;
    
%-------------------------------------------------

i=1; fn{i}= filename;
eval(['load ',filename,';']); load(fn{i}) ;  ip=find(fn{i} == '.')-1; c{i}=eval(fn{i}(1:ip));
nC=size(c{1},2);
nR=size(c{1});
if nC == 8
    % --- old version without CQ
    c{i}(:,nC+1)=  c{i}(:,2)-c{i}(:,nC);  % Cx-Fore
    cqv=0;
else
    % new version, with CQ
    c{i}(:,nC+1)= c{i}(:,nC);
    c{i}(:,nC)=  c{i}(:,2)-c{i}(:,nC-1);  % Cx-Fore
    cqv=1;
    cq=c{i}(:,nC+1);
end
jCX_FORE=nC+1-cqv;jCX_BASE=jCX_FORE-1;
alow=min(c{i}(:,1)); ahigh=max(c{i}(:,1)); axp=[alow-5,ahigh+5];
% --------------plot -----------
%
mr(1)={'-ob'}; mrl(1)={'-b'};
%
% convegrnce rating. index is QC+2
%
conr{2}='or'; conr{4}='om';conr{6}='oc'; conr{8}='ob';
conr{12}='og'; conr{1}='ok';

figure(1);
%
  subplot(311)

    kvar=2;
  plot(c{i}(:,1),c{i}(:,kvar),mrl{i},axp,[0,0],'-k','LineWidth',2,'MarkerSize',6)
  h1=gca; set(h1,'FontSize',14);
  title(['Configuration: ',Confcase]);
  grid on; xlabel('\alpha [^o]'); ylabel('CX');
  if cqv==1,polaraPlot_setQualMap;end


%
  subplot(312)

  plot(c{i}(:,1),c{i}(:,jCX_FORE),mr{i},axp,[0,0],'-k','LineWidth',2,'MarkerSize',6);
        
  h1=gca; set(h1,'FontSize',14);
  grid on; xlabel('\alpha [^o]'); ylabel('CX-FORE');
  %
%
  subplot(313)

  plot(c{i}(:,1),c{i}(:,jCX_BASE),mr{i},axp,[0,0],'-k','LineWidth',2,'MarkerSize',6);
        
  h1=gca; set(h1,'FontSize',14);
  grid on; xlabel('\alpha [^o]'); ylabel('CX-Base');
 
figure(2);
%
%----------- Normal plane -------------
%
%

  subplot(211)
    kvar=3;
  plot(c{i}(:,1),c{i}(:,3),mrl{i},axp,[0,0],'-k','LineWidth',2,'MarkerSize',6)
  h1=gca; set(h1,'FontSize',14);
  title(['Configuration: ',Confcase]);
  grid on; xlabel('\alpha [^o]'); ylabel('-CZ = CNOR');
  if cqv==1,polaraPlot_setQualMap;end


%
  subplot(212)

  plot(c{i}(:,1),-c{i}(:,7),mr{i},axp,[0,0],'-k','LineWidth',2,'MarkerSize',6);
        
  h1=gca; set(h1,'FontSize',14);
  grid on; xlabel('\alpha [^o]'); ylabel('Cm');
  %
%
figure(3);
%
%--- yaw and roll
%
  subplot(311)

    kvar=-4;
  plot(c{i}(:,1),-c{i}(:,4),mrl{i},axp,[0,0],'-k','LineWidth',2,'MarkerSize',6)
  h1=gca; set(h1,'FontSize',14);
  title(['Configuration: ',Confcase]);
  grid on; xlabel('\alpha [^o]'); ylabel('CY (side-force)');
  if cqv==1,polaraPlot_setQualMap;end


%
  subplot(312)

  plot(c{i}(:,1),-c{i}(:,6),mr{i},axp,[0,0],'-k','LineWidth',2,'MarkerSize',6);
        
  h1=gca; set(h1,'FontSize',14);
  grid on; xlabel('\alpha [^o]'); ylabel('Cn (yaw)');
  %
%
  subplot(313)

  plot(c{i}(:,1),-c{i}(:,5),mr{i},axp,[0,0],'-k','LineWidth',2,'MarkerSize',6);
        
  h1=gca; set(h1,'FontSize',14);
 
  grid on; xlabel('\alpha [^o]'); ylabel('Cl (roll)');
 
