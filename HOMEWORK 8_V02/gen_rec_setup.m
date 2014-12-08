% This script set parameter and performs intial calculations for
% waveform level modulation of a Permanent magnet AC machine

% Written by Maryam Alibeik
% Department of Electrical and Computer Engineering
% Purdue School of Engineering and Technology, IUPUI
% 723 West Michigan St. SL 113
% Indianapolis, IN 46202
% PH: 317-370-5153
% Email: malibeik@purdue.edu

P.Ll=1.19e-03;                          %DC-Link Inductor Inductance (H)
P.rl=0.32;                              %DC-Link series resistance (Ohm)
P.rs=382e-03;                           %Leakage Inductor series resistance in the synchronous machine (Ohm)
P.Lls=1.12e-03;                         %Leakage Inductance (H)
P.Lmq=24.9e-03;                         %Magnetizing q-axis Inductance (H)
P.Lmd=39.3e-03;                         %Magnetizing D-axis Inductance (H)
P.rd1=140;                              %D-axis leakage resistance (Ohm)                              
P.rd2=1.19e+03;                         %D-axis leakage resistance (Ohm) 
P.rd3=1.58;                             %D-axis leakage resistance (Ohm)
P.rfd=112e-03;                          %field series resistance (Ohm)
P.Lld1=9.87e-03;                        %D-axis leakage inductance (H)
P.Lld2=4.91e-03;                        %D-axis leakage inductance (H)
P.Lld3=4.52e-03;                        %D-axis leakage inductance (H)
P.Llfd=1.53e-03;                        %leakage field inductance    
P.rq1=5.07;                             %Q-axis leakage resistance (Ohm)
P.rq2=1.06;                             %Q-axis leakage resistance (Ohm)
P.rq3=447e-03;                          %Q-axis leakage resistance (Ohm)
P.pole=4;                               %number of poles
P.Llq1=4.21e-03;                        %Q-axis leakage inductance (H)
P.Llq2=3.5e-03;                         %Q-axis leakage inductance (H)
P.Llq3=26.2e-03;                        %Q-axis leakage inductance (H)