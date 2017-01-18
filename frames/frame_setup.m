function [analysis,adjoint,PhiPhiTr,no_of_frames]=frame_setup(z,par,Nsc)
% Sets up variables for main optimization algorithm
% [analysis,adjoint,PhiPhiT,no_of_frames]=frame_setup(z,par)
% analysis/adjoint - handle to analysis/synthesi function
% z - degraded image
% par - structure containing parameters (ie. No.of iteration, Type of frame)

no_of_frames=length(par.frame);
if ~exist('Nsc','var'), Nsc = 3; end % number of scales
analysis = cell(no_of_frames,1);
adjoint = cell(no_of_frames,1);
PhiPhiT = cell(no_of_frames,1);
for m = 1:no_of_frames
    switch par.frame{m}
        case 'TV4' % Total Variation - zda se, ze je tam nejaky problem?
            analysis{m}=@(z)analysisTV4(z);
            adjoint{m}=@(pyr)adjointTV4(pyr);
            PhiPhiT{m} = ppsf2otf([0 -1 0;-1 4 -1;0 -1 0],size(z));
        case 'TV4i' % Total Variation - zda se, ze je tam nejaky problem?
            analysis{m}=@(z)analysisTV4i(z);
            adjoint{m}=@(pyr)adjointTV4i(pyr);
            PhiPhiT{m} = ppsf2otf([0 -1 0;-1 4 -1;0 -1 0],size(z));
        case 'TV6' % Total Variation - zda se, ze je tam nejaky problem?
            analysis{m}=@(z)analysisTV6(z);
            adjoint{m}=@(pyr)adjointTV6(pyr);
            lapl=zeros(3,3,3);
            lapl(1:3,2,2)=[-1 2 -1];
            lapl=lapl+shiftdim(lapl,1)+shiftdim(lapl,2);
            PhiPhiT{m} = single(ppsf2otf(lapl,size(z)));        % POZOR !!! je single kvuli pameti
        case 'TV8' % - ZATIM MOC NEFUNGUJE
            analysis{m}=@(z)analysisTV8(z);
            adjoint{m}=@(pyr)adjointTV8(pyr);
            kernel = conv2([-1 1;0 0],flipc([-1 1;0 0]),'full')+...
                conv2([-1 0; 1 0],flipc([-1 0; 1 0]),'full')+...
                conv2([-1 0;0 1],flipc([-1 0;0 1]),'full')+...
                conv2([0 -1;1 0],flipc([-1 0;0 1]),'full'); %\Phi\Phi^T
            PhiPhiT{m} = ppsf2otf(kernel,size(z));
        case 'DT-CWT' % Dual Tree Complex Wavelet Transform
            Nsc=5;
            par2 = 'near_sym_a'; par3 = 'qshift_a'; % Pyramid parameters
            analysis{m}=@(z)analysisDTCWT(z,Nsc,par2,par3);  % Simoncelli
            adjoint{m}=@(pyr)synthesisDTCWT(pyr,par2,par3); % Simoncelli
            PhiPhiT{m} = 1;
        case 'TIHP' % Translation Invariant Haar Pyramid
            Nsc = 3;
            analysis{m}=@(z)analysisTIHP(z,Nsc); % Portilla
            adjoint{m}=@(pyr)synthesisTIHP(pyr); % Portilla
            PhiPhiT{m} = 1;
        case 'mdwt' % classical dwt
            h = daubcqf(2,'mid');
            warning('Wrong now');
            analysis{m}=@(z)mdwt(z,h,Nsc);
            adjoint{m}=@(pyr)midwt(pyr,h,Nsc);
            PhiPhiT{m} = 1;
        case 'rdwt' % redundant dwt
            h = daubcqf(2,'mid');
            analysis{m}=@(z)mrdwt(z,h,Nsc);
            adjoint{m}=@(pyr)mirdwt(pyr,h,Nsc);
            PhiPhiT{m} = 1;
        case 'TIHP2' % redundant dwt - analougous to Javier's implementation - stale funguje z nejakeho duvodu hur
            Nsc = 7;
            h = daubcqf(2,'mid');
            analysis{m}=@(z)mrdwt2(z,h,Nsc);
            adjoint{m}=@(pyr)mirdwt2(pyr,h,Nsc);
            PhiPhiT{m} = 1; %?? mozna neni korektni - zalezi jestli je to normalizovany frame
        case 'Haar' % Discrete Wavelet Transform - Matlab Wavelet Toolbox
            Nsc=5;
            par2 = 'db1'; % Pyramid parameters, type of wavelet, see help of wfilters.m
            analysis{m}=@(z)analysisDWT(z,Nsc,par2);
            adjoint{m}=@(pyr)synthesisDWT(pyr,par2);
            PhiPhiT{m} = 1;
        case 'DB2' % Discrete Wavelet Transform - Matlab Wavelet Toolbox
            Nsc=5;
            par2 = 'db2'; % Pyramid parameters, type of wavelet, see help of wfilters.m
            analysis{m}=@(z)analysisDWT(z,Nsc,par2);
            adjoint{m}=@(pyr)synthesisDWT(pyr,par2);
            PhiPhiT{m} = 1;
        case 'DB3' % Discrete Wavelet Transform - Matlab Wavelet Toolbox
            Nsc=5;
            par2 = 'db3'; % Pyramid parameters, type of wavelet, see help of wfilters.m
            analysis{m}=@(z)analysisDWT(z,Nsc,par2);
            adjoint{m}=@(pyr)synthesisDWT(pyr,par2);
            PhiPhiT{m} = 1;
        case 'LearnedD' % Learned dictionary
            analysis{m}=@(z)analysisLearnedD(z,par.dictionary);
            adjoint{m}=@(pyr)synthesisLearnedD(pyr,par.dictionary);
            PhiPhiT{m} = 1;
        otherwise
            error('Unknown frame type');
    end
end
PhiPhiTr = PhiPhiT{1};
for i = 2:no_of_frames
    PhiPhiTr = PhiPhiTr + PhiPhiT{m};
end
PhiPhiTr = PhiPhiTr / no_of_frames; % sum_i \Phi_i\Phi^T_i
