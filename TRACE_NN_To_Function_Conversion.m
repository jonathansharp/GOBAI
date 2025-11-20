% Network to Function converter.    
VerboseTF=1;
for Property=7:7;
    if     Property==1; VName='Preformed_O'; Equations=1; FN='ESPER_PP_Preformed_O.mat';
    elseif Property==2; VName='Preformed_P'; Equations=1; FN='ESPER_PP_Preformed_P.mat';
    elseif Property==3; VName='Preformed_Si'; Equations=1; FN='ESPER_PP_Preformed_Si.mat';
    elseif Property==4; VName='Preformed_N'; Equations=1; FN='ESPER_PP_Preformed_N.mat';
    elseif Property==5; VName='Preformed_TA'; Equations=1; FN='ESPER_PP_Preformed_TA.mat';
    elseif Property==6; VName='EstT_Temperature'; Equations=1; FN='ESPER_EstT_Temperature.mat';
    elseif Property==7; VName='SFs'; Equations=1; FN='../MonteCarlo/NNets/TRACE_NN_IG_Sandborn_Test_0_SF.mat';
    end
    % Loading the files needed
    load(FN);
    
    for Equation=1:1 % stepping through all 16 equations
        for n=1:4 % A committee of 4 neural networks is used.
            % Separate neural networks are used for the Arctic/Atlantic and
            % the rest of the ocean (more info below).
            cc(Nets.(strcat('Eqn',num2str(Equation))).Atl.(strcat('Net',num2str(n))),horzcat('./TRACE_NN_Functions/TRACETest_',VName,'_',num2str(Equation),'_Atl_',num2str(n)));
            genFunction(Nets.(strcat('Eqn',num2str(Equation))).Other.(strcat('Net',num2str(n))),horzcat('./TRACE_NN_Functions/TRACETest_',VName,'_',num2str(Equation),'_Other_',num2str(n)));
        end
        if VerboseTF==true; % Updating on progress if desired
            disp(horzcat(num2str(floor(toc)),' seconds elapsed... now using eqn. ',num2str(Equation)))
        end
    end
end