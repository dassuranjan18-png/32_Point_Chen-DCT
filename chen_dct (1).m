function Feature1 = chen_dct(x)
format long e;

numRuns = 37388;                   % Number of runs

[N, M] = size(x);        % N = 32, M = number of frames
Feature1 = zeros(N, M, 'single');

for runIdx = 1:numRuns
f = single(x(:, runIdx));





% Quantize inputs once by DCB (as if converter present at inputs)
for n = 1:32
    f(n) = dcb(f(n));
end

% ---------- First Stage FDCT (fully unrolled with DCB on every operand) ----------
F_stage1 = zeros(32, numRuns, 'single');

F_stage1(0+1,runIdx)  = dcb( dcb(f(0+1))  + dcb(f(31+1))  );
F_stage1(31+1,runIdx) = dcb( dcb(f(0+1))  - dcb(f(31+1))  );
F_stage1(1+1,runIdx)  = dcb( dcb(f(1+1))  + dcb(f(30+1))  );
F_stage1(30+1,runIdx) = dcb( dcb(f(1+1))  - dcb(f(30+1))  );
F_stage1(2+1,runIdx)  = dcb( dcb(f(2+1))  + dcb(f(29+1))  );
F_stage1(29+1,runIdx) = dcb( dcb(f(2+1))  - dcb(f(29+1))  );
F_stage1(3+1,runIdx)  = dcb( dcb(f(3+1))  + dcb(f(28+1))  );
F_stage1(28+1,runIdx) = dcb( dcb(f(3+1))  - dcb(f(28+1))  );
F_stage1(4+1,runIdx)  = dcb( dcb(f(4+1))  + dcb(f(27+1)) );
F_stage1(27+1,runIdx) = dcb( dcb(f(4+1))  - dcb(f(27+1)) );
F_stage1(5+1,runIdx)  = dcb( dcb(f(5+1))  + dcb(f(26+1)) );
F_stage1(26+1,runIdx) = dcb( dcb(f(5+1))  - dcb(f(26+1)) );
F_stage1(6+1,runIdx)  = dcb( dcb(f(6+1))  + dcb(f(25+1)) );
F_stage1(25+1,runIdx) = dcb( dcb(f(6+1))  - dcb(f(25+1)) );
F_stage1(7+1,runIdx)  = dcb( dcb(f(7+1))  + dcb(f(24+1)) );
F_stage1(24+1,runIdx) = dcb( dcb(f(7+1))  - dcb(f(24+1)) );
F_stage1(8+1,runIdx)  = dcb( dcb(f(8+1))  + dcb(f(23+1)) );
F_stage1(23+1,runIdx) = dcb( dcb(f(8+1))  - dcb(f(23+1)) );
F_stage1(9+1,runIdx)  = dcb( dcb(f(9+1))  + dcb(f(22+1)) );
F_stage1(22+1,runIdx) = dcb( dcb(f(9+1))  - dcb(f(22+1)) );
F_stage1(10+1,runIdx) = dcb( dcb(f(10+1)) + dcb(f(21+1)) );
F_stage1(21+1,runIdx) = dcb( dcb(f(10+1)) - dcb(f(21+1)) );
F_stage1(11+1,runIdx) = dcb( dcb(f(11+1)) + dcb(f(20+1)) );
F_stage1(20+1,runIdx) = dcb( dcb(f(11+1)) - dcb(f(20+1)) );
F_stage1(12+1,runIdx) = dcb( dcb(f(12+1)) + dcb(f(19+1)) );
F_stage1(19+1,runIdx) = dcb( dcb(f(12+1)) - dcb(f(19+1)) );
F_stage1(13+1,runIdx) = dcb( dcb(f(13+1)) + dcb(f(18+1)) );
F_stage1(18+1,runIdx) = dcb( dcb(f(13+1)) - dcb(f(18+1)) );
F_stage1(14+1,runIdx) = dcb( dcb(f(14+1)) + dcb(f(17+1)) );
F_stage1(17+1,runIdx) = dcb( dcb(f(14+1)) - dcb(f(17+1)) );
F_stage1(15+1,runIdx) = dcb( dcb(f(15+1)) + dcb(f(16+1)) );
F_stage1(16+1,runIdx) = dcb( dcb(f(15+1)) - dcb(f(16+1)) );




% ---------- Second Stage FDCT (fully unrolled with DCB) ----------
% Using C = cos(pi/4) = S = sin(pi/4) = sqrt(2)/2
F_stage2 = zeros(32, numRuns, 'single');
C  = single(cos(pi/4));   % +0.70710677
Cn = single(-cos(pi/4));  % -0.70710677



% --- Upper branch (indices 0..15): pair (0,15), (1,14), ..., (7,8)
% ===== Stage 2 =====
F_stage2(0+1,runIdx)  = dcb( dcb(F_stage1(0+1,runIdx))  + dcb(F_stage1(15+1,runIdx)) );
F_stage2(15+1,runIdx) = dcb( dcb(F_stage1(0+1,runIdx))  - dcb(F_stage1(15+1,runIdx)) );
F_stage2(1+1,runIdx)  = dcb( dcb(F_stage1(1+1,runIdx))  + dcb(F_stage1(14+1,runIdx)) );
F_stage2(14+1,runIdx) = dcb( dcb(F_stage1(1+1,runIdx))  - dcb(F_stage1(14+1,runIdx)) );
F_stage2(2+1,runIdx)  = dcb( dcb(F_stage1(2+1,runIdx))  + dcb(F_stage1(13+1,runIdx)) );
F_stage2(13+1,runIdx) = dcb( dcb(F_stage1(2+1,runIdx))  - dcb(F_stage1(13+1,runIdx)) );
F_stage2(3+1,runIdx)  = dcb( dcb(F_stage1(3+1,runIdx))  + dcb(F_stage1(12+1,runIdx)) );
F_stage2(12+1,runIdx) = dcb( dcb(F_stage1(3+1,runIdx))  - dcb(F_stage1(12+1,runIdx)) );
F_stage2(4+1,runIdx)  = dcb( dcb(F_stage1(4+1,runIdx))  + dcb(F_stage1(11+1,runIdx)) );
F_stage2(11+1,runIdx) = dcb( dcb(F_stage1(4+1,runIdx))  - dcb(F_stage1(11+1,runIdx)) );
F_stage2(5+1,runIdx)  = dcb( dcb(F_stage1(5+1,runIdx))  + dcb(F_stage1(10+1,runIdx)) );
F_stage2(10+1,runIdx) = dcb( dcb(F_stage1(5+1,runIdx))  - dcb(F_stage1(10+1,runIdx)) );
F_stage2(6+1,runIdx)  = dcb( dcb(F_stage1(6+1,runIdx))  + dcb(F_stage1(9+1,runIdx)) );
F_stage2(9+1,runIdx)  = dcb( dcb(F_stage1(6+1,runIdx))  - dcb(F_stage1(9+1,runIdx)) );
F_stage2(7+1,runIdx)  = dcb( dcb(F_stage1(7+1,runIdx))  + dcb(F_stage1(8+1,runIdx)) );
F_stage2(8+1,runIdx)  = dcb( dcb(F_stage1(7+1,runIdx))  - dcb(F_stage1(8+1,runIdx)) );

F_stage2(16+1,runIdx) = F_stage1(16+1,runIdx);
F_stage2(31+1,runIdx) = F_stage1(31+1,runIdx);
F_stage2(17+1,runIdx) = F_stage1(17+1,runIdx);
F_stage2(30+1,runIdx) = F_stage1(30+1,runIdx);
F_stage2(18+1,runIdx) = F_stage1(18+1,runIdx);
F_stage2(29+1,runIdx) = F_stage1(29+1,runIdx);
F_stage2(19+1,runIdx) = F_stage1(19+1,runIdx);
F_stage2(28+1,runIdx) = F_stage1(28+1,runIdx);

F_stage2(20+1,runIdx) = dcb( dcb(dcb(F_stage1(20+1,runIdx))* dcb(Cn)) + dcb(dcb(F_stage1(27+1,runIdx))* dcb(C)) );
F_stage2(27+1,runIdx) = dcb( dcb(dcb(F_stage1(20+1,runIdx))* dcb(C)) + dcb(dcb(F_stage1(27+1,runIdx))* dcb(C)) );
F_stage2(21+1,runIdx) = dcb( dcb(dcb(F_stage1(21+1,runIdx))* dcb(Cn)) + dcb(dcb(F_stage1(26+1,runIdx))* dcb(C)) );
F_stage2(26+1,runIdx) = dcb( dcb(dcb(F_stage1(21+1,runIdx))* dcb(C)) + dcb(dcb(F_stage1(26+1,runIdx))* dcb(C)) );
F_stage2(22+1,runIdx) = dcb( dcb(dcb(F_stage1(22+1,runIdx))* dcb(Cn)) + dcb(dcb(F_stage1(25+1,runIdx))* dcb(C)) );
F_stage2(25+1,runIdx) = dcb( dcb(dcb(F_stage1(22+1,runIdx))* dcb(C)) + dcb(dcb(F_stage1(25+1,runIdx))* dcb(C)) );
F_stage2(23+1,runIdx) = dcb( dcb(dcb(F_stage1(23+1,runIdx))* dcb(Cn)) + dcb(dcb(F_stage1(24+1,runIdx))* dcb(C)) );
F_stage2(24+1,runIdx) = dcb( dcb(dcb(F_stage1(23+1,runIdx))* dcb(C)) + dcb(dcb(F_stage1(24+1,runIdx))* dcb(C)) );








% ---------- 3rd Stage FDCT (fully unrolled with DCB) ----------
% Using C = cos(pi/4) = S = sin(pi/4) = sqrt(2)/2
F_stage3 = zeros(32, numRuns, 'single');


% --- Upper half (indices 0..15) : pairs (0,7), (1,6), (2,5), (3,4)
% ===== Stage 3 =====
F_stage3(0+1,runIdx) = dcb( dcb(F_stage2(0+1,runIdx)) + dcb(F_stage2(7+1,runIdx)) );
F_stage3(7+1,runIdx) = dcb( dcb(F_stage2(0+1,runIdx)) - dcb(F_stage2(7+1,runIdx)) );
F_stage3(1+1,runIdx) = dcb( dcb(F_stage2(1+1,runIdx)) + dcb(F_stage2(6+1,runIdx)) );
F_stage3(6+1,runIdx) = dcb( dcb(F_stage2(1+1,runIdx)) - dcb(F_stage2(6+1,runIdx)) );
F_stage3(2+1,runIdx) = dcb( dcb(F_stage2(2+1,runIdx)) + dcb(F_stage2(5+1,runIdx)) );
F_stage3(5+1,runIdx) = dcb( dcb(F_stage2(2+1,runIdx)) - dcb(F_stage2(5+1,runIdx)) );
F_stage3(3+1,runIdx) = dcb( dcb(F_stage2(3+1,runIdx)) + dcb(F_stage2(4+1,runIdx)) );
F_stage3(4+1,runIdx) = dcb( dcb(F_stage2(3+1,runIdx)) - dcb(F_stage2(4+1,runIdx)) );

F_stage3(8+1,runIdx)  = F_stage2(8+1,runIdx);
F_stage3(9+1,runIdx)  = F_stage2(9+1,runIdx);
F_stage3(14+1,runIdx) = F_stage2(14+1,runIdx);
F_stage3(15+1,runIdx) = F_stage2(15+1,runIdx);

F_stage3(10+1,runIdx) = dcb( dcb(dcb(F_stage2(10+1,runIdx))*dcb(Cn)) + dcb(dcb(F_stage2(13+1,runIdx))*dcb(C)) );
F_stage3(13+1,runIdx) = dcb( dcb(dcb(F_stage2(10+1,runIdx))*dcb(C))  + dcb(dcb(F_stage2(13+1,runIdx))*dcb(C)) );
F_stage3(11+1,runIdx) = dcb( dcb(dcb(F_stage2(11+1,runIdx))*dcb(Cn)) + dcb(dcb(F_stage2(12+1,runIdx))*dcb(C)) );
F_stage3(12+1,runIdx) = dcb( dcb(dcb(F_stage2(11+1,runIdx))*dcb(C))  + dcb(dcb(F_stage2(12+1,runIdx))*dcb(C)) );

F_stage3(16+1,runIdx) = dcb( dcb(F_stage2(16+1,runIdx)) + dcb(F_stage2(23+1,runIdx)) );
F_stage3(23+1,runIdx) = dcb( dcb(F_stage2(16+1,runIdx)) - dcb(F_stage2(23+1,runIdx)) );
F_stage3(17+1,runIdx) = dcb( dcb(F_stage2(17+1,runIdx)) + dcb(F_stage2(22+1,runIdx)) );
F_stage3(22+1,runIdx) = dcb( dcb(F_stage2(17+1,runIdx)) - dcb(F_stage2(22+1,runIdx)) );
F_stage3(18+1,runIdx) = dcb( dcb(F_stage2(18+1,runIdx)) + dcb(F_stage2(21+1,runIdx)) );
F_stage3(21+1,runIdx) = dcb( dcb(F_stage2(18+1,runIdx)) - dcb(F_stage2(21+1,runIdx)) );
F_stage3(19+1,runIdx) = dcb( dcb(F_stage2(19+1,runIdx)) + dcb(F_stage2(20+1,runIdx)) );
F_stage3(20+1,runIdx) = dcb( dcb(F_stage2(19+1,runIdx)) - dcb(F_stage2(20+1,runIdx)) );

F_stage3(24+1,runIdx) = dcb( - dcb(F_stage2(24+1,runIdx)) + dcb(F_stage2(31+1,runIdx)) );
F_stage3(31+1,runIdx) = dcb( dcb(F_stage2(24+1,runIdx))  + dcb(F_stage2(31+1,runIdx)) );
F_stage3(25+1,runIdx) = dcb( - dcb(F_stage2(25+1,runIdx)) + dcb(F_stage2(30+1,runIdx)) );
F_stage3(30+1,runIdx) = dcb( dcb(F_stage2(25+1,runIdx))  + dcb(F_stage2(30+1,runIdx)) );
F_stage3(26+1,runIdx) = dcb( - dcb(F_stage2(26+1,runIdx)) + dcb(F_stage2(29+1,runIdx)) );
F_stage3(29+1,runIdx) = dcb( dcb(F_stage2(26+1,runIdx))  + dcb(F_stage2(29+1,runIdx)) );
F_stage3(27+1,runIdx) = dcb( - dcb(F_stage2(27+1,runIdx)) + dcb(F_stage2(28+1,runIdx)) );
F_stage3(28+1,runIdx) = dcb( dcb(F_stage2(27+1,runIdx))  + dcb(F_stage2(28+1,runIdx)) );






% ---------- Fourth Stage FDCT (fully unrolled with DCB) ----------
% Using C = cos(pi/4) = S = sin(pi/4) = sqrt(2)/2
F_stage4 = zeros(32, numRuns, 'single');
C8   = single(cos(pi/8));     % ≈ 0.9238795
S8   = single(sin(pi/8));     % ≈ 0.38268343
C3_8 = single(cos(3*pi/8));   % ≈ 0.38268343
S3_8 = single(sin(3*pi/8));   % ≈ 0.9238795
C8n   = single(-C8);
S8n   = single(-S8);
C3_8n = single(-C3_8);
S3_8n = single(-S3_8);

F_stage4(0+1, runIdx)  = dcb( dcb(F_stage3(0+1,runIdx)) + dcb(F_stage3(3+1,runIdx)) );
F_stage4(3+1, runIdx)  = dcb( dcb(F_stage3(0+1,runIdx)) - dcb(F_stage3(3+1,runIdx)) );
F_stage4(1+1, runIdx)  = dcb( dcb(F_stage3(1+1,runIdx)) + dcb(F_stage3(2+1,runIdx)) );
F_stage4(2+1, runIdx)  = dcb( dcb(F_stage3(1+1,runIdx)) - dcb(F_stage3(2+1,runIdx)) );

F_stage4(8+1, runIdx)  = dcb( dcb(F_stage3(8+1,runIdx)) + dcb(F_stage3(11+1,runIdx)) );
F_stage4(11+1, runIdx) = dcb( dcb(F_stage3(8+1,runIdx)) - dcb(F_stage3(11+1,runIdx)) );
F_stage4(9+1, runIdx)  = dcb( dcb(F_stage3(9+1,runIdx)) + dcb(F_stage3(10+1,runIdx)) );
F_stage4(10+1, runIdx) = dcb( dcb(F_stage3(9+1,runIdx)) - dcb(F_stage3(10+1,runIdx)) );

F_stage4(12+1, runIdx) = dcb(- dcb(F_stage3(12+1,runIdx)) + dcb(F_stage3(15+1,runIdx)) );
F_stage4(15+1, runIdx) = dcb( dcb(F_stage3(12+1,runIdx)) + dcb(F_stage3(15+1,runIdx)) );
F_stage4(13+1, runIdx) = dcb(- dcb(F_stage3(13+1,runIdx)) + dcb(F_stage3(14+1,runIdx)) );
F_stage4(14+1, runIdx) = dcb( dcb(F_stage3(13+1,runIdx)) + dcb(F_stage3(14+1,runIdx)) );

F_stage4(5+1, runIdx)  = dcb( dcb(dcb(F_stage3(5+1,runIdx))*dcb(Cn)) + dcb(dcb(F_stage3(6+1,runIdx))*dcb(C)) );
F_stage4(6+1, runIdx)  = dcb( dcb(dcb(F_stage3(5+1,runIdx))*dcb(C))  + dcb(dcb(F_stage3(6+1,runIdx))*dcb(C)) );

F_stage4(18+1, runIdx) = dcb( dcb(F_stage3(18+1,runIdx))* dcb(C8n) + dcb(F_stage3(29+1,runIdx))* dcb(S8)) ;
F_stage4(29+1, runIdx) = dcb( dcb(dcb(F_stage3(18+1,runIdx))* dcb(C3_8)) + dcb(dcb(F_stage3(29+1,runIdx))* dcb(S3_8)));
F_stage4(19+1, runIdx) = dcb( dcb(dcb(F_stage3(19+1,runIdx))* dcb(C8n)) + dcb(dcb(F_stage3(28+1,runIdx))* dcb(S8))) ;
F_stage4(28+1, runIdx) = dcb( dcb(dcb(F_stage3(19+1,runIdx))* dcb(C3_8)) + dcb(dcb(F_stage3(28+1,runIdx))* dcb(S3_8)) );
F_stage4(20+1, runIdx) = dcb( dcb(dcb(F_stage3(20+1,runIdx))* dcb(S8n)) + dcb(dcb(F_stage3(27+1,runIdx))* dcb(C8n))) ;
F_stage4(27+1, runIdx) = dcb( dcb(dcb(F_stage3(20+1,runIdx))* dcb(S3_8n)) + dcb(dcb(F_stage3(27+1,runIdx))* dcb(C3_8)) );
F_stage4(21+1, runIdx) = dcb( dcb(dcb(F_stage3(21+1,runIdx))* dcb(S8n)) + dcb(dcb(F_stage3(26+1,runIdx))* dcb(C8n))) ;
F_stage4(26+1, runIdx) = dcb( dcb(dcb(F_stage3(21+1,runIdx))* dcb(S3_8n)) + dcb(dcb(F_stage3(26+1,runIdx))* dcb(C3_8)) );

F_stage4(4+1, runIdx)  = F_stage3(4+1,runIdx);
F_stage4(7+1, runIdx)  = F_stage3(7+1,runIdx);
F_stage4(16+1, runIdx) = F_stage3(16+1,runIdx);
F_stage4(17+1, runIdx) = F_stage3(17+1,runIdx);
F_stage4(30+1, runIdx) = F_stage3(30+1,runIdx);
F_stage4(31+1, runIdx) = F_stage3(31+1,runIdx);
F_stage4(22+1, runIdx) = F_stage3(22+1,runIdx);
F_stage4(23+1, runIdx) = F_stage3(23+1,runIdx);
F_stage4(24+1, runIdx) = F_stage3(24+1,runIdx);
F_stage4(25+1, runIdx) = F_stage3(25+1,runIdx);




% ---------- 5th Stage FDCT (fully unrolled with DCB) ----------
% Using C = cos(pi/4) = S = sin(pi/4) = sqrt(2)/2
F_stage5 = zeros(32, numRuns, 'single');


F_stage5(16+1, runIdx) = dcb( dcb(F_stage4(16+1,runIdx)) + dcb(F_stage4(19+1,runIdx)) );
F_stage5(19+1, runIdx) = dcb( dcb(F_stage4(16+1,runIdx)) - dcb(F_stage4(19+1,runIdx)) );
F_stage5(17+1, runIdx) = dcb( dcb(F_stage4(17+1,runIdx)) + dcb(F_stage4(18+1,runIdx)) );
F_stage5(18+1, runIdx) = dcb( dcb(F_stage4(17+1,runIdx)) - dcb(F_stage4(18+1,runIdx)) );

F_stage5(20+1, runIdx) = dcb( -dcb(F_stage4(20+1,runIdx)) + dcb(F_stage4(23+1,runIdx)) );
F_stage5(23+1, runIdx) = dcb( dcb(F_stage4(20+1,runIdx)) + dcb(F_stage4(23+1,runIdx)) );
F_stage5(21+1, runIdx) = dcb(- dcb(F_stage4(21+1,runIdx)) + dcb(F_stage4(22+1,runIdx)) );
F_stage5(22+1, runIdx) = dcb( dcb(F_stage4(21+1,runIdx)) + dcb(F_stage4(22+1,runIdx)) );

F_stage5(24+1, runIdx) = dcb(- dcb(F_stage4(24+1,runIdx)) + dcb(F_stage4(27+1,runIdx)) );
F_stage5(27+1, runIdx) = dcb( dcb(F_stage4(24+1,runIdx)) + dcb(F_stage4(27+1,runIdx)) );
F_stage5(25+1, runIdx) = dcb(- dcb(F_stage4(25+1,runIdx)) + dcb(F_stage4(25+1,runIdx)) );
F_stage5(26+1, runIdx) = dcb( dcb(F_stage4(25+1,runIdx)) + dcb(F_stage4(25+1,runIdx)) );

F_stage5(28+1, runIdx) = dcb( dcb(F_stage4(28+1,runIdx)) + dcb(F_stage4(31+1,runIdx)) );
F_stage5(31+1, runIdx) = dcb( dcb(F_stage4(28+1,runIdx)) - dcb(F_stage4(31+1,runIdx)) );
F_stage5(29+1, runIdx) = dcb( dcb(F_stage4(29+1,runIdx)) + dcb(F_stage4(30+1,runIdx)) );
F_stage5(30+1, runIdx) = dcb( dcb(F_stage4(29+1,runIdx)) - dcb(F_stage4(30+1,runIdx)) );

F_stage5(4+1, runIdx)  = dcb( dcb(F_stage4(4+1,runIdx)) + dcb(F_stage4(5+1,runIdx)) );
F_stage5(5+1, runIdx)  = dcb( dcb(F_stage4(4+1,runIdx)) - dcb(F_stage4(5+1,runIdx)) );
F_stage5(6+1, runIdx)  = dcb( dcb(F_stage4(6+1,runIdx)) + dcb(F_stage4(7+1,runIdx)) );
F_stage5(7+1, runIdx)  = dcb( dcb(F_stage4(6+1,runIdx)) - dcb(F_stage4(7+1,runIdx)) );

F_stage5(0+1, runIdx)  = dcb( dcb(dcb(F_stage4(0+1,runIdx))*dcb(C)) + dcb(dcb(F_stage4(1+1,runIdx))*dcb(C)) );
F_stage5(1+1, runIdx)  = dcb( dcb(dcb(F_stage4(0+1,runIdx))*dcb(C))  + dcb(dcb(F_stage4(1+1,runIdx))*dcb(Cn)) );
F_stage5(2+1, runIdx)  = dcb( dcb(dcb(F_stage4(2+1,runIdx))* dcb(S8)) + dcb(dcb(F_stage4(3+1,runIdx))* dcb(C8))) ;
F_stage5(3+1, runIdx)  = dcb( dcb(dcb(F_stage4(2+1,runIdx))* dcb(S3_8n)) + dcb(dcb(F_stage4(3+1,runIdx))* dcb(C3_8)) );

F_stage5(9+1, runIdx)  = dcb( dcb(dcb(F_stage4(9+1,runIdx))* dcb(C8n)) + dcb(dcb(F_stage4(14+1,runIdx))* dcb(S8))) ;
F_stage5(14+1, runIdx) = dcb( dcb(dcb(F_stage4(9+1,runIdx))* dcb(C3_8)) + dcb(dcb(F_stage4(14+1,runIdx))* dcb(S3_8) ));
F_stage5(10+1, runIdx) = dcb( dcb(dcb(F_stage4(10+1,runIdx))* dcb(S8n)) + dcb(dcb(F_stage4(27+1,runIdx))* dcb(C8n))) ;
F_stage5(13+1, runIdx) = dcb( dcb(dcb(F_stage4(10+1,runIdx))* dcb(S3_8)) + dcb(dcb(F_stage4(27+1,runIdx))* dcb(C3_8)) );

F_stage5(8+1, runIdx)  = F_stage4(8+1,runIdx);
F_stage5(15+1, runIdx) = F_stage4(15+1,runIdx);
F_stage5(11+1, runIdx) = F_stage4(11+1,runIdx);
F_stage5(12+1, runIdx) = F_stage4(12+1,runIdx);






% ---------- 6th Stage FDCT (fully unrolled with DCB) ----------
% Using C = cos(pi/4) = S = sin(pi/4) = sqrt(2)/2
F_stage6 = zeros(32, numRuns, 'single');

C16 = single(cos(pi/16));
S16 = single(sin(pi/16));
C3_16 = single(cos(3*pi/16));
S3_16 = single(sin(3*pi/16));
C5_16 = single(cos(5*pi/16));
S5_16 = single(sin(5*pi/16));
C7_16 = single(cos(7*pi/16));
S7_16 = single(sin(7*pi/16));
C16n   = single(-cos(pi/16));
S16n   = single(-sin(pi/16));
C3_16n = single(-cos(3*pi/16));
S3_16n = single(-sin(3*pi/16));
C5_16n = single(-cos(5*pi/16));
S5_16n = single(-sin(5*pi/16));
C7_16n = single(-cos(7*pi/16));
S7_16n = single(-sin(7*pi/16));

F_stage6(4+1, runIdx)  = dcb( dcb(dcb(F_stage5(4+1,runIdx))*dcb(S16)) + dcb(dcb(F_stage5(7+1,runIdx))*dcb(C16)) );
F_stage6(7+1, runIdx)  = dcb( dcb(dcb(F_stage5(4+1,runIdx))*dcb(S7_16n))  + dcb(dcb(F_stage5(7+1,runIdx))*dcb(C7_16)) );
F_stage6(5+1, runIdx)  = dcb( dcb(dcb(F_stage5(5+1,runIdx))* dcb(S5_16)) + dcb(dcb(F_stage5(6+1,runIdx))* dcb(C5_16))) ;
F_stage6(6+1, runIdx)  = dcb( dcb(dcb(F_stage5(5+1,runIdx))* dcb(S3_16n)) + dcb(dcb(F_stage5(6+1,runIdx))* dcb(C3_16)) );
F_stage6(8+1, runIdx)  = dcb( dcb(F_stage5(8+1,runIdx)) + dcb(F_stage5(9+1,runIdx)) );
F_stage6(9+1, runIdx)  = dcb( dcb(F_stage5(8+1,runIdx)) - dcb(F_stage5(9+1,runIdx)) );
F_stage6(10+1, runIdx) = dcb( -dcb(F_stage5(10+1,runIdx)) + dcb(F_stage5(11+1,runIdx)) );
F_stage6(11+1, runIdx) = dcb( dcb(F_stage5(10+1,runIdx)) + dcb(F_stage5(11+1,runIdx)) );
F_stage6(12+1, runIdx) = dcb( dcb(F_stage5(12+1,runIdx)) + dcb(F_stage5(13+1,runIdx)) );
F_stage6(13+1, runIdx) = dcb( dcb(F_stage5(12+1,runIdx)) - dcb(F_stage5(13+1,runIdx)) );
F_stage6(14+1, runIdx) = dcb( -dcb(F_stage5(14+1,runIdx)) + dcb(F_stage5(15+1,runIdx)) );
F_stage6(15+1, runIdx) = dcb( dcb(F_stage5(14+1,runIdx)) + dcb(F_stage5(15+1,runIdx)) );

F_stage6(17+1, runIdx) = dcb( dcb(dcb(F_stage5(17+1,runIdx))*dcb(C16n)) + dcb(dcb(F_stage5(30+1,runIdx))*dcb(S16)) );
F_stage6(30+1, runIdx) = dcb( dcb(dcb(F_stage5(17+1,runIdx))*dcb(C7_16))  + dcb(dcb(F_stage5(30+1,runIdx))*dcb(S7_16)) );
F_stage6(18+1, runIdx) = dcb( dcb(dcb(F_stage5(18+1,runIdx))* dcb(S16n)) + dcb(dcb(F_stage5(29+1,runIdx))* dcb(C16n))) ;
F_stage6(29+1, runIdx) = dcb( dcb(dcb(F_stage5(18+1,runIdx))* dcb(S7_16n)) + dcb(dcb(F_stage5(29+1,runIdx))* dcb(C7_16)) );

F_stage6(21+1, runIdx) = dcb( dcb(dcb(F_stage5(21+1,runIdx))* dcb(C5_16n)) + dcb(dcb(F_stage5(26+1,runIdx))* dcb(S5_16))) ;
F_stage6(26+1, runIdx) = dcb( dcb(dcb(F_stage5(21+1,runIdx))* dcb(C3_16)) + dcb(dcb(F_stage5(26+1,runIdx))* dcb(S3_16)) );
F_stage6(22+1, runIdx) = dcb( dcb(dcb(F_stage5(22+1,runIdx))* dcb(S5_16n)) + dcb(dcb(F_stage5(25+1,runIdx))* dcb(C5_16n))) ;
F_stage6(25+1, runIdx) = dcb( dcb(dcb(F_stage5(22+1,runIdx))* dcb(S3_16n)) + dcb(dcb(F_stage5(25+1,runIdx))* dcb(C3_16)) );

F_stage6(16+1, runIdx) = F_stage5(16+1,runIdx);
F_stage6(19+1, runIdx) = F_stage5(19+1,runIdx);
F_stage6(20+1, runIdx) = F_stage5(20+1,runIdx);
F_stage6(23+1, runIdx) = F_stage5(23+1,runIdx);
F_stage6(24+1, runIdx) = F_stage5(24+1,runIdx);
F_stage6(27+1, runIdx) = F_stage5(27+1,runIdx);
F_stage6(28+1, runIdx) = F_stage5(28+1,runIdx);
F_stage6(29+1, runIdx) = F_stage5(29+1,runIdx);
F_stage6(30+1, runIdx) = F_stage5(30+1,runIdx);
F_stage6(31+1, runIdx) = F_stage5(31+1,runIdx);












% ---------- 7th Stage FDCT (fully unrolled with DCB) ----------
% Using C = cos(pi/4) = S = sin(pi/4) = sqrt(2)/2
F_stage7 = zeros(32, numRuns, 'single');


C32   =  single( cos(1*pi/32));
C32n  = -single( cos(1*pi/32));
S32   =  single( sin(1*pi/32));
S32n  = -single( sin(1*pi/32));
C3_32   =  single( cos(3*pi/32));
C3_32n  = -single( cos(3*pi/32));
S3_32   =  single( sin(3*pi/32));
S3_32n  = -single( sin(3*pi/32));
C5_32   =  single( cos(5*pi/32));
C5_32n  = -single( cos(5*pi/32));
S5_32   =  single( sin(5*pi/32));
S5_32n  = -single( sin(5*pi/32));
C7_32   =  single( cos(7*pi/32));
C7_32n  = -single( cos(7*pi/32));
S7_32   =  single( sin(7*pi/32));
S7_32n  = -single( sin(7*pi/32));
C9_32   =  single( cos(9*pi/32));
C9_32n  = -single( cos(9*pi/32));
S9_32   =  single( sin(9*pi/32));
S9_32n  = -single( sin(9*pi/32));
C11_32  =  single( cos(11*pi/32));
C11_32n = -single( cos(11*pi/32));
S11_32  =  single( sin(11*pi/32));
S11_32n = -single( sin(11*pi/32));
C13_32  =  single( cos(13*pi/32));
C13_32n = -single( cos(13*pi/32));
S13_32  =  single( sin(13*pi/32));
S13_32n = -single( sin(13*pi/32));
C15_32  =  single( cos(15*pi/32));
C15_32n = -single( cos(15*pi/32));
S15_32  =  single( sin(15*pi/32));
S15_32n = -single( sin(15*pi/32));



F_stage7(16+1, runIdx) = dcb( dcb(F_stage6(16+1,runIdx)) + dcb(F_stage6(17+1,runIdx)) );
F_stage7(17+1, runIdx) = dcb( dcb(F_stage6(16+1,runIdx)) - dcb(F_stage6(17+1,runIdx)) );
F_stage7(18+1, runIdx) = dcb( -dcb(F_stage6(18+1,runIdx)) + dcb(F_stage6(19+1,runIdx)) );
F_stage7(19+1, runIdx) = dcb( dcb(F_stage6(18+1,runIdx)) + dcb(F_stage6(19+1,runIdx)) );
F_stage7(20+1, runIdx) = dcb( dcb(F_stage6(20+1,runIdx)) + dcb(F_stage6(21+1,runIdx)) );
F_stage7(21+1, runIdx) = dcb( dcb(F_stage6(20+1,runIdx)) - dcb(F_stage6(21+1,runIdx)) );
F_stage7(22+1, runIdx) = dcb( -dcb(F_stage6(22+1,runIdx)) + dcb(F_stage6(23+1,runIdx)) );
F_stage7(23+1, runIdx) = dcb( dcb(F_stage6(22+1,runIdx)) + dcb(F_stage6(23+1,runIdx)) );
F_stage7(24+1, runIdx) = dcb( dcb(F_stage6(24+1,runIdx)) + dcb(F_stage6(25+1,runIdx)) );
F_stage7(25+1, runIdx) = dcb( dcb(F_stage6(24+1,runIdx)) - dcb(F_stage6(25+1,runIdx)) );
F_stage7(26+1, runIdx) = dcb( -dcb(F_stage6(26+1,runIdx)) + dcb(F_stage6(27+1,runIdx)) );
F_stage7(27+1, runIdx) = dcb( dcb(F_stage6(26+1,runIdx)) + dcb(F_stage6(27+1,runIdx)) );
F_stage7(28+1, runIdx) = dcb( dcb(F_stage6(28+1,runIdx)) + dcb(F_stage6(29+1,runIdx)) );
F_stage7(29+1, runIdx) = dcb( dcb(F_stage6(28+1,runIdx)) - dcb(F_stage6(29+1,runIdx)) );
F_stage7(30+1, runIdx) = dcb( -dcb(F_stage6(30+1,runIdx)) + dcb(F_stage6(31+1,runIdx)) );
F_stage7(31+1, runIdx) = dcb( dcb(F_stage6(30+1,runIdx)) + dcb(F_stage6(31+1,runIdx)) );

F_stage7(8+1, runIdx)  = dcb( dcb(dcb(F_stage6(8+1,runIdx))*dcb(S32)) + dcb(dcb(F_stage6(15+1,runIdx))*dcb(C3_32)) );
F_stage7(15+1, runIdx) = dcb( dcb(dcb(F_stage6(8+1,runIdx))*dcb(S15_32n))  + dcb(dcb(F_stage6(15+1,runIdx))*dcb(C15_32)) );
F_stage7(9+1, runIdx)  = dcb( dcb(dcb(F_stage6(9+1,runIdx))* dcb(S9_32)) + dcb(dcb(F_stage6(14+1,runIdx))* dcb(C9_32))) ;
F_stage7(14+1, runIdx) = dcb( dcb(dcb(F_stage6(9+1,runIdx))* dcb(S7_32n)) + dcb(dcb(F_stage6(14+1,runIdx))* dcb(C7_32)) );
F_stage7(10+1, runIdx) = dcb( dcb(dcb(F_stage6(10+1,runIdx))* dcb(S5_32)) + dcb(dcb(F_stage6(13+1,runIdx))* dcb(C5_32))) ;
F_stage7(13+1, runIdx) = dcb( dcb(dcb(F_stage6(10+1,runIdx))* dcb(C11_32n)) + dcb(dcb(F_stage6(13+1,runIdx))* dcb(C11_32)) );
F_stage7(11+1, runIdx) = dcb( dcb(dcb(F_stage6(11+1,runIdx))* dcb(S13_32)) + dcb(dcb(F_stage6(12+1,runIdx))* dcb(C13_32))) ;
F_stage7(12+1, runIdx) = dcb( dcb(dcb(F_stage6(11+1,runIdx))* dcb(S3_32n)) + dcb(dcb(F_stage6(12+1,runIdx))* dcb(C3_32)));













% ---------- 8th Stage FDCT (fully unrolled with DCB) ----------
% Using C = cos(pi/4) = S = sin(pi/4) = sqrt(2)/2
F_stage8 = zeros(32, numRuns, 'single');

C64  =  single( cos(pi/64));   S64  =  single( sin(pi/64));   C64n  = single(-cos(pi/64));   S64n  = single(-sin(pi/64));
C3_64  =  single( cos(3*pi/64)); S3_64  =  single( sin(3*pi/64)); C3_64n  = single(-cos(3*pi/64)); S3_64n  = single(-sin(3*pi/64));
C5_64  =  single( cos(5*pi/64)); S5_64  =  single( sin(5*pi/64)); C5_64n  = single(-cos(5*pi/64)); S5_64n  = single(-sin(5*pi/64));
C7_64  =  single( cos(7*pi/64)); S7_64  =  single( sin(7*pi/64)); C7_64n  = single(-cos(7*pi/64)); S7_64n  = single(-sin(7*pi/64));
C9_64  =  single( cos(9*pi/64)); S9_64  =  single( sin(9*pi/64)); C9_64n  = single(-cos(9*pi/64)); S9_64n  = single(-sin(9*pi/64));
C11_64 =  single(cos(11*pi/64)); S11_64 =  single(sin(11*pi/64)); C11_64n = single(-cos(11*pi/64));S11_64n = single(-sin(11*pi/64));
C13_64 =  single(cos(13*pi/64)); S13_64 =  single(sin(13*pi/64)); C13_64n = single(-cos(13*pi/64));S13_64n = single(-sin(13*pi/64));
C15_64 =  single(cos(15*pi/64)); S15_64 =  single(sin(15*pi/64)); C15_64n = single(-cos(15*pi/64));S15_64n = single(-sin(15*pi/64));
C17_64 =  single(cos(17*pi/64)); S17_64 =  single(sin(17*pi/64)); C17_64n = single(-cos(17*pi/64));S17_64n = single(-sin(17*pi/64));
C19_64 =  single(cos(19*pi/64)); S19_64 =  single(sin(19*pi/64)); C19_64n = single(-cos(19*pi/64));S19_64n = single(-sin(19*pi/64));
C21_64 =  single(cos(21*pi/64)); S21_64 =  single(sin(21*pi/64)); C21_64n = single(-cos(21*pi/64));S21_64n = single(-sin(21*pi/64));
C23_64 =  single(cos(23*pi/64)); S23_64 =  single(sin(23*pi/64)); C23_64n = single(-cos(23*pi/64));S23_64n = single(-sin(23*pi/64));
C25_64 =  single(cos(25*pi/64)); S25_64 =  single(sin(25*pi/64)); C25_64n = single(-cos(25*pi/64));S25_64n = single(-sin(25*pi/64));
C27_64 =  single(cos(27*pi/64)); S27_64 =  single(sin(27*pi/64)); C27_64n = single(-cos(27*pi/64));S27_64n = single(-sin(27*pi/64));
C29_64 =  single(cos(29*pi/64)); S29_64 =  single(sin(29*pi/64)); C29_64n = single(-cos(29*pi/64));S29_64n = single(-sin(29*pi/64));
C31_64 =  single(cos(31*pi/64)); S31_64 =  single(sin(31*pi/64)); C31_64n = single(-cos(31*pi/64));S31_64n = single(-sin(31*pi/64));



F_stage8(16+1, runIdx) = dcb( dcb(dcb(F_stage7(16+1,runIdx))*dcb(S64))     + dcb(dcb(F_stage7(31+1,runIdx))*dcb(C64)) );
F_stage8(31+1, runIdx) = dcb( dcb(dcb(F_stage7(16+1,runIdx))*dcb(S31_64n)) + dcb(dcb(F_stage7(31+1,runIdx))*dcb(C31_64)) );
F_stage8(17+1, runIdx) = dcb( dcb(dcb(F_stage7(17+1,runIdx))* dcb(S17_64))  + dcb(dcb(F_stage7(30+1,runIdx))* dcb(C17_64))) ;
F_stage8(30+1, runIdx) = dcb( dcb(dcb(F_stage7(17+1,runIdx))* dcb(S15_64n)) + dcb(dcb(F_stage7(30+1,runIdx))* dcb(C15_64)) );
F_stage8(18+1, runIdx) = dcb( dcb(dcb(F_stage7(18+1,runIdx))* dcb(S9_64)) + dcb(dcb(F_stage7(29+1,runIdx))* dcb(C9_64))) ;
F_stage8(29+1, runIdx) = dcb( dcb(dcb(F_stage7(18+1,runIdx))* dcb(S23_64n)) + dcb(dcb(F_stage7(29+1,runIdx))* dcb(C23_64)));
F_stage8(19+1, runIdx) = dcb( dcb(dcb(F_stage7(19+1,runIdx))* dcb(S25_64)) + dcb(dcb(F_stage7(28+1,runIdx))* dcb(C25_64))) ;
F_stage8(28+1, runIdx) = dcb( dcb(dcb(F_stage7(19+1,runIdx))* dcb(S7_64n)) + dcb(dcb(F_stage7(28+1,runIdx))* dcb(C7_64)));
F_stage8(20+1, runIdx) = dcb( dcb(dcb(F_stage7(20+1,runIdx))*dcb(S5_64))     + dcb(dcb(F_stage7(27+1,runIdx))*dcb(C5_64)) );
F_stage8(27+1, runIdx) = dcb( dcb(dcb(F_stage7(20+1,runIdx))*dcb(S27_64n)) + dcb(dcb(F_stage7(27+1,runIdx))*dcb(C27_64)) );
F_stage8(21+1, runIdx) = dcb( dcb(dcb(F_stage7(21+1,runIdx))* dcb(S21_64))  + dcb(dcb(F_stage7(26+1,runIdx))* dcb(C21_64))) ;
F_stage8(26+1, runIdx) = dcb( dcb(dcb(F_stage7(21+1,runIdx))* dcb(S11_64n)) + dcb(dcb(F_stage7(26+1,runIdx))* dcb(C11_64)));
F_stage8(22+1, runIdx) = dcb( dcb(dcb(F_stage7(22+1,runIdx))* dcb(S13_64)) + dcb(dcb(F_stage7(25+1,runIdx))* dcb(C13_64))) ;
F_stage8(25+1, runIdx) = dcb( dcb(dcb(F_stage7(22+1,runIdx))* dcb(S19_64n)) + dcb(dcb(F_stage7(25+1,runIdx))* dcb(C19_64)));
F_stage8(23+1, runIdx) = dcb( dcb(dcb(F_stage7(23+1,runIdx))* dcb(S29_64)) + dcb(dcb(F_stage7(24+1,runIdx))* dcb(C29_64))) ;
F_stage8(24+1, runIdx) = dcb( dcb(dcb(F_stage7(23+1,runIdx))* dcb(S3_64n)) + dcb(dcb(F_stage7(24+1,runIdx))* dcb(C3_64)) );















% ---------- Final output FDCT (fully unrolled with DCB) ----------

dctII = zeros(32, numRuns, 'single');




dctII(0+1, runIdx)  = F_stage5(0+1, runIdx)  * 0.1767767;
dctII(1+1, runIdx)  = F_stage8(16+1, runIdx) * 0.25;
dctII(2+1, runIdx)  = F_stage7(8+1, runIdx)  * 0.25;
dctII(3+1, runIdx)  = F_stage8(24+1, runIdx) * 0.25;
dctII(4+1, runIdx)  = F_stage6(4+1, runIdx)  * 0.25;
dctII(5+1, runIdx)  = F_stage8(20+1, runIdx) * 0.25;
dctII(6+1, runIdx)  = F_stage7(12+1, runIdx) * 0.25;
dctII(7+1, runIdx)  = F_stage8(28+1, runIdx) * 0.25;
dctII(8+1, runIdx)  = F_stage5(2+1, runIdx)  * 0.25;
dctII(9+1, runIdx)  = F_stage8(18+1, runIdx) * 0.25;
dctII(10+1, runIdx) = F_stage7(10+1, runIdx) * 0.25;
dctII(11+1, runIdx) = F_stage8(26+1, runIdx) * 0.25;
dctII(12+1, runIdx) = F_stage6(6+1, runIdx)  * 0.25;
dctII(13+1, runIdx) = F_stage8(22+1, runIdx) * 0.25;
dctII(14+1, runIdx) = F_stage7(14+1, runIdx) * 0.25;
dctII(15+1, runIdx) = F_stage8(30+1, runIdx) * 0.25;
dctII(16+1, runIdx) = F_stage5(1+1, runIdx)  * 0.25;
dctII(17+1, runIdx) = F_stage8(17+1, runIdx) * 0.25;
dctII(18+1, runIdx) = F_stage7(9+1, runIdx)  * 0.25;
dctII(19+1, runIdx) = F_stage8(25+1, runIdx) * 0.25;
dctII(20+1, runIdx) = F_stage6(5+1, runIdx)  * 0.25;
dctII(21+1, runIdx) = F_stage8(21+1, runIdx) * 0.25;
dctII(22+1, runIdx) = F_stage7(13+1, runIdx) * 0.25;
dctII(23+1, runIdx) = F_stage8(29+1, runIdx) * 0.25;
dctII(24+1, runIdx) = F_stage5(3+1, runIdx)  * 0.25;
dctII(25+1, runIdx) = F_stage8(19+1, runIdx) * 0.25;
dctII(26+1, runIdx) = F_stage7(11+1, runIdx) * 0.25;
dctII(27+1, runIdx) = F_stage8(27+1, runIdx) * 0.25;
dctII(28+1, runIdx) = F_stage6(7+1, runIdx)  * 0.25;
dctII(29+1, runIdx) = F_stage8(23+1, runIdx) * 0.25;
dctII(30+1, runIdx) = F_stage7(15+1, runIdx) * 0.25;
dctII(31+1, runIdx) = F_stage8(31+1, runIdx) * 0.25;








Feature1(:, runIdx) = dctII(:, runIdx); 
end

save('SD.mat', 'Feature1');
% ===========================
% === Local Functions Below =
% ===========================

%function y = dcb(x,runIdx)
    % DCB bypass mode: just return input as single precision
    %y = single(x);
%end



function y = dcb(x)
     %DCB: convert single -> CFP(18-bit) -> single (introduces quantization error)
    x = double(x);
    y = arrayfun(@(v) unpack16(pack16(v)), x, 'UniformOutput', true);
    y = single(y);
    end
    

function bits = pack18(v)
    % Pack a double scalar into uint32 containing the 18-bit CFP pattern:
    % [bit17 sign][bits16..10 exponent (7 bits)][bits9..0 mantissa (10 bits)]
    % Bias = 63. Denormal support. Round-to-nearest mantissa.
    v = double(v);
    if isnan(v), bits = uint32(0); return; end
    s = uint32(v < 0);
    av = abs(v);
    if av == 0
        bits = bitor(bitshift(s,17), uint32(0));
        return;
    end
    bias = 63;
    e = floor(log2(av));           % unbiased exponent
    normv = av / 2.^e;             % in [1,2)
    E_stored = e + bias;
    if E_stored <= 0
        % Denormal range: exponent field = 0
        mant_raw = round(av ./ 2.^(1 - bias) .* 2.^10);
        mant_raw = min(max(mant_raw, 0), 1023);
        E_field = uint32(0);
        M_field = uint32(mant_raw);
    else
        % Normalized
        mant = normv - 1;
        mant_raw = round(mant .* 2.^10);
        if mant_raw >= 1024
            mant_raw = 0;
            E_stored = E_stored + 1;
        end
        if E_stored >= 127
            E_field = uint32(127);
            M_field = uint32(1023);
        else
            E_field = uint32(E_stored);
            M_field = uint32(min(max(mant_raw,0),1023));
        end
    end
    bits = bitor(bitshift(s,17), bitor(bitshift(E_field,10), M_field));
end

function v = unpack18(bits)
    % Unpack uint32 containing 18-bit CFP pattern back to double scalar
    b = uint32(bits);
    s = bitshift(b, -17);
    E_field = bitand(bitshift(b, -10), uint32(127));
    M_field = bitand(b, uint32(1023));
    bias = 63;
    if E_field == 0
        mant = double(M_field) / 2.^10;
        vmag = mant * 2.^(1 - bias);
    else
        mant = 1 + double(M_field) / 2.^10;
        vmag = mant * 2.^(double(E_field) - bias);
    end
    v = ((-1) .^ double(s)) .* vmag;
end




function bits = pack16(v)
    % Pack a double scalar into uint16 containing 16-bit CEP pattern:
    % [bit15 sign][bits14..8 exponent (7 bits)][bits7..0 mantissa (8 bits)]
    % Bias = 63. Denormal support. Round-to-nearest mantissa.

    v = double(v);
    if isnan(v), bits = uint16(0); return; end
    s = uint16(v < 0);
    av = abs(v);

    if av == 0
        bits = bitor(bitshift(s,15), uint16(0));
        return;
    end

    bias = 63;
    e = floor(log2(av));          % unbiased exponent
    normv = av / 2.^e;            % in [1,2)
    E_stored = e + bias;

    if E_stored <= 0
        % Denormal
        mant_raw = round(av ./ 2.^(1 - bias) .* 2.^8);
        mant_raw = min(max(mant_raw, 0), 255);
        E_field = uint16(0);
        M_field = uint16(mant_raw);
    else
        % Normalized
        mant = normv - 1;
        mant_raw = round(mant .* 2.^8); % 8-bit mantissa
        if mant_raw >= 256
            mant_raw = 0;
            E_stored = E_stored + 1;
        end
        if E_stored >= 127
            % Overflow → Inf (max exponent & mantissa)
            E_field = uint16(127);
            M_field = uint16(255);
        else
            E_field = uint16(E_stored);
            M_field = uint16(min(max(mant_raw,0),255));
        end
    end

    bits = bitor(bitshift(s,15), bitor(bitshift(E_field,8), M_field));
end

function v = unpack16(bits)
    % Unpack uint16 containing 16-bit CEP pattern back to double scalar
    b = uint16(bits);
    s = bitshift(b, -15);
    E_field = bitand(bitshift(b, -8), uint16(127));
    M_field = bitand(b, uint16(255));
    bias = 63;

    if E_field == 0
        % Denormal
        mant = double(M_field) / 2.^8;
        vmag = mant * 2.^(1 - bias);
    else
        mant = 1 + double(M_field) / 2.^8;
        vmag = mant * 2.^(double(E_field) - bias);
    end

    v = ((-1) .^ double(s)) .* vmag;
end




function bits = pack14(v)
    % Pack a double scalar into uint16 containing the 14-bit CEP pattern:
    % [bit13 sign][bits12..6 exponent (7 bits)][bits5..0 mantissa (6 bits)]
    % Bias = 63. Denormal support. Round-to-nearest mantissa.
    v = double(v);
    if isnan(v), bits = uint16(0); return; end
    
    s = uint16(v < 0);
    av = abs(v);
    if av == 0
        bits = bitshift(s,13); % just the sign in bit 13
        return;
    end

    bias = 63;
    e = floor(log2(av));        % unbiased exponent
    normv = av / 2.^e;          % in [1,2)
    E_stored = e + bias;
    
    if E_stored <= 0
        % Denormal: exponent field = 0
        mant_raw = round(av ./ 2.^(1 - bias) .* 2.^6);
        mant_raw = min(max(mant_raw, 0), 63);
        E_field = uint16(0);
        M_field = uint16(mant_raw);
    else
        % Normalized
        mant = normv - 1;
        mant_raw = round(mant .* 2.^6);
        if mant_raw >= 64
            mant_raw = 0;
            E_stored = E_stored + 1;
        end
        if E_stored >= 127
            % Overflow -> max exponent and mantissa
            E_field = uint16(127);
            M_field = uint16(63);
        else
            E_field = uint16(E_stored);
            M_field = uint16(min(max(mant_raw,0),63));
        end
    end
    
    bits = bitor(bitshift(s,13), bitor(bitshift(E_field,6), M_field));
end

function v = unpack14(bits)
    % Unpack uint16 containing 14-bit CEP pattern back to double scalar
    b = uint16(bits);
    s = bitshift(b, -13);
    E_field = bitand(bitshift(b, -6), uint16(127));
    M_field = bitand(b, uint16(63));
    bias = 63;
    
    if E_field == 0
        % Denormal
        mant = double(M_field) / 2.^6;
        vmag = mant * 2.^(1 - bias);
    else
        mant = 1 + double(M_field) / 2.^6;
        vmag = mant * 2.^(double(E_field) - bias);
    end
    
    v = ((-1) .^ double(s)) .* vmag;
end
end