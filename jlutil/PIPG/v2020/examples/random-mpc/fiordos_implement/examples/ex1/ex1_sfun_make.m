function ex1_sfun_make(varargin)
%     floatprec   - precision of float data type
%                   allowed values: 'double', 'single'
%                   default: 'double'

% 
% This file is generated by FiOrdOs, a program licensed under GPL
% by copyright holder Automatic Control Laboratory, ETH Zurich.
% 
% If you are interested in using this file commercially,
% please contact the copyright holder.
% 

ip=inputParser();
ip.addParamValue('floatprec','double', @(val)(ismember(val,{'double','single'})));
ip.parse(varargin{:})
use_realtype_single = strcmp(ip.Results.floatprec,'single');

if use_realtype_single
    mex -v -DUSE_REALTYPE_SINGLE ex1_sfun.c
else
    mex -v ex1_sfun.c
end


if bdIsLoaded('ex1_sfun_lib'),
    bdclose('ex1_sfun_lib');
end
if exist('ex1_sfun_lib.mdl','file')==4
    delete('ex1_sfun_lib.mdl');
end
if exist('ex1_sfun_lib.slx','file')==4
    delete('ex1_sfun_lib.slx');
end
new_system('ex1_sfun_lib','Library');

add_block('built-in/SubSystem','ex1_sfun_lib/ex1_sfun');
set_param('ex1_sfun_lib/ex1_sfun','MaskPromptString','algo.maxit source:| algo.maxit value [1x1]:|algo.init source:| algo.init value [10x1]:|algo.stopgEps source:| algo.stopgEps value [1x1]:|algo.stopgStride [1x1]:|Sample time (>0):');
set_param('ex1_sfun_lib/ex1_sfun','MaskVariables','d_algo_maxit=@1;s_algo_maxit=@2;d_algo_init=@3;s_algo_init=@4;d_algo_stopgEps=@5;s_algo_stopgEps=@6;s_algo_stopgStride=@7;ts=@8;');
set_param('ex1_sfun_lib/ex1_sfun','MaskStyleString','popup(internal|external),edit,popup(internal|external),edit,popup(internal|external),edit,edit,edit');
set_param('ex1_sfun_lib/ex1_sfun','MaskValueString','internal|5000|internal|[0;0;0;0;0;0;0;0;0;0]|internal|0.001|1|1');
set_param('ex1_sfun_lib/ex1_sfun','MaskTunableValueString','off,off,off,off,off,off,off,off');
set_param('ex1_sfun_lib/ex1_sfun','MaskVisibilityString','on,on,on,on,on,on,on,on');
set_param('ex1_sfun_lib/ex1_sfun','MaskTabNameString','general,general,init,init,algo.stop,algo.stop,algo.stop,general');
set_param('ex1_sfun_lib/ex1_sfun','MaskCallBackString','ex1_sfun_lib_aux(''cb'',gcb,1)||ex1_sfun_lib_aux(''cb'',gcb,3)||ex1_sfun_lib_aux(''cb'',gcb,5)|||');

add_block('built-in/S-Function','ex1_sfun_lib/ex1_sfun/ex1_sfun0');
set_param('ex1_sfun_lib/ex1_sfun/ex1_sfun0','Position',[250 50 550 170]);
set_param('ex1_sfun_lib/ex1_sfun','Position',[50 50 350 170]);
set_param('ex1_sfun_lib/ex1_sfun/ex1_sfun0','Parameters','s_algo_stopgStride,ts');
set_param('ex1_sfun_lib/ex1_sfun/ex1_sfun0','FunctionName','ex1_sfun');
drawCmds={
    'port_label(''input'',1,''g [10x1]'')'
    'port_label(''input'',2,''c [1x1]'')'
    'port_label(''input'',3,''\italgo.init [10x1]'',''texmode'',''on'')'
    'port_label(''input'',4,''\italgo.maxit [1x1]'',''texmode'',''on'')'
    'port_label(''input'',5,''\italgo.stopgEps [1x1]'',''texmode'',''on'')'
    'port_label(''output'',1,''x [10x1]'')'
    'port_label(''output'',2,''f [1x1]'')'
    'port_label(''output'',3,''iter [1x1]'')'
    'port_label(''output'',4,''exitflag [1x1]'')'
    'disp(''FiOrdOs'')'
};
set_param('ex1_sfun_lib/ex1_sfun/ex1_sfun0','MaskDisplay',char(drawCmds));

add_block('built-in/Inport','ex1_sfun_lib/ex1_sfun/g [10x1]','Position',[50 50 80 64]);
add_line('ex1_sfun_lib/ex1_sfun','g [10x1]/1','ex1_sfun0/1');
add_block('built-in/Inport','ex1_sfun_lib/ex1_sfun/c [1x1]','Position',[50 89 80 103]);
add_line('ex1_sfun_lib/ex1_sfun','c [1x1]/1','ex1_sfun0/2');

add_block('built-in/Constant','ex1_sfun_lib/ex1_sfun/algo.init [10x1]','Position',[50 128 80 142],'Value','s_algo_init');
add_line('ex1_sfun_lib/ex1_sfun','algo.init [10x1]/1','ex1_sfun0/3');
add_block('built-in/Constant','ex1_sfun_lib/ex1_sfun/algo.maxit [1x1]','Position',[50 167 80 181],'Value','s_algo_maxit');
add_line('ex1_sfun_lib/ex1_sfun','algo.maxit [1x1]/1','ex1_sfun0/4');
add_block('built-in/Constant','ex1_sfun_lib/ex1_sfun/algo.stopgEps [1x1]','Position',[50 206 80 220],'Value','s_algo_stopgEps');
add_line('ex1_sfun_lib/ex1_sfun','algo.stopgEps [1x1]/1','ex1_sfun0/5');

add_block('built-in/Outport','ex1_sfun_lib/ex1_sfun/x [10x1]','Position',[700 50 730 64]);
add_line('ex1_sfun_lib/ex1_sfun','ex1_sfun0/1','x [10x1]/1');
add_block('built-in/Outport','ex1_sfun_lib/ex1_sfun/f [1x1]','Position',[700 89 730 103]);
add_line('ex1_sfun_lib/ex1_sfun','ex1_sfun0/2','f [1x1]/1');
add_block('built-in/Outport','ex1_sfun_lib/ex1_sfun/iter [1x1]','Position',[700 128 730 142]);
add_line('ex1_sfun_lib/ex1_sfun','ex1_sfun0/3','iter [1x1]/1');
add_block('built-in/Outport','ex1_sfun_lib/ex1_sfun/exitflag [1x1]','Position',[700 167 730 181]);
add_line('ex1_sfun_lib/ex1_sfun','ex1_sfun0/4','exitflag [1x1]/1');

set_param('ex1_sfun_lib/ex1_sfun','MaskSelfModifiable','on');
set_param('ex1_sfun_lib/ex1_sfun','MaskInitialization','ex1_sfun_lib_aux(''init'',gcb)');

save_system('ex1_sfun_lib','ex1_sfun_lib');
set_param('ex1_sfun_lib', 'Lock','on');
open_system('ex1_sfun_lib');

end