% Simple script to tabulate case data into a structure.

% Define the strucutre table.
T=struct('Name',[], 'Filename',[], 'Mc',[], 'Mk',[], 'Ts',NaT(0),'Te',NaT(0),'md',[]);
i=1;


% Basel.
T(i).Name='Basel';
T(i).Filename='data/Basel.mat';
T(i).Mc=+0.9;
T(i).dMb=0.1;
T(i).dMd=0.0;
T(i).Ts( 1)=datetime(2006,12,02, 18,18,33); T(i).Te( 1)=datetime(2006,12,08, 11,22,45); T(i).md( 1)=0;
T(i).Tf( 1)=datetime(2006,12,20, 01,01,01);
T(i).dtd=0;
T(i).S=[1];
i=i+1;

% SSFS-1993.
T(i).Name='SSFS93';
T(i).Filename='/Users/rschultz/Desktop/papers/old/Trailing/IS-bath/data/processed/SSFS1993.mat';
T(i).Mc=-1.5;
T(i).dMb=0.1;
T(i).dMd=0.00;
T(i).Ts( 1)=datetime(1993,09,01, 20,30,38); T(i).Te( 1)=datetime(1993,10,16, 13,34,28); T(i).md( 1)=0;
T(i).Tf( 1)=datetime(1993,09,22, 20,41,01);
T(i).dtd=0;
T(i).S=[1];
i=i+1;

% SSFS-2000.
T(i).Name='SSFS00';
T(i).Filename='/Users/rschultz/Desktop/papers/old/Trailing/IS-bath/data/processed/SSFS2000.mat';
T(i).Mc=-0.5;
T(i).dMb=0.2;
T(i).dMd=0.2;
T(i).Ts( 1)=datetime(2000,06,30, 18,39,59); T(i).Te( 1)=datetime(2000,07,06, 16,40,10); T(i).md( 1)=0;
T(i).Tf( 1)=datetime(2000,07,11, 05,55,19);
T(i).dtd=1;
T(i).S=[1];
i=i+1;

% SSFS-2003.
T(i).Name='SSFS03';
T(i).Filename='/Users/rschultz/Desktop/papers/old/Trailing/IS-bath/data/processed/SSFS2003.mat';
T(i).Mc=+0.1;
T(i).dMb=0.2;
T(i).dMd=0.2;
T(i).Ts( 1)=datetime(2003,05,27, 16,15,06); T(i).Te( 1)=datetime(2003,07,16, 19,56,56); T(i).md( 1)=0;
T(i).Tf( 1)=datetime(2003,06,11, 08,06,00);
T(i).dtd=1;
T(i).S=[1];
i=i+1;

% SSFS-2004.
T(i).Name='SSFS04';
T(i).Filename='/Users/rschultz/Desktop/papers/old/Trailing/IS-bath/data/processed/SSFS2004.mat';
T(i).Mc=-0.8;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts( 1)=datetime(2004,09,13, 08,03,56); T(i).Te( 1)=datetime(2004,09,16, 20,00,00); T(i).md( 1)=0;
T(i).Tf( 1)=datetime(2004,09,24, 12,18,36);
T(i).dtd=1;
T(i).S=[1];
i=i+1;

% SSFS-2005.
T(i).Name='SSFS05';
T(i).Filename='/Users/rschultz/Desktop/papers/old/Trailing/IS-bath/data/processed/SSFS2005.mat';
T(i).Mc=-0.2;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts( 1)=datetime(2005,02,07, 18,00,32); T(i).Te( 1)=datetime(2005,02,11, 18,20,42); T(i).md( 1)=0;
T(i).Ts( 2)=datetime(2005,02,22, 15,29,28); T(i).Te( 2)=datetime(2005,02,25, 15,40,48); T(i).md( 2)=0;
T(i).Ts( 3)=datetime(2005,03,02, 12,42,06); T(i).Te( 3)=datetime(2005,03,05, 06,13,13); T(i).md( 3)=0;
T(i).Ts( 4)=datetime(2005,03,13, 20,39,18); T(i).Te( 4)=datetime(2005,03,16, 20,55,42); T(i).md( 4)=0;
T(i).Tf( 1)=datetime(2005,02,16, 09,00,00);
T(i).dtd=0;
T(i).S=[1];
i=i+1;

% Cooper Basin - HAB4.
T(i).Name='CB4';
T(i).Filename='/Users/rschultz/Desktop/papers/old/Trailing/IS-bath/data/processed/CB_HAB4.mat';
T(i).Mc=+0.9;
T(i).dMb=0.3;
T(i).dMd=0.1;
T(i).Ts( 1)=datetime(2012,11,13, 09,25,55); T(i).Te( 1)=datetime(2012,11,30, 12,59,42); T(i).md( 1)=0;
T(i).Tf( 1)=datetime(2099,01,01, 01,00,00);
T(i).dtd=1;
T(i).S=[1];
i=i+1;

% Helsinki St1 2018.
T(i).Name='St1-2018';
T(i).Filename='/Users/rschultz/Desktop/papers/old/Trailing/IS-bath/data/processed/St1-2018.mat';
T(i).Mc=0.0;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(2018,06,04, 06,26,43); T(i).Te(1)=datetime(2018,06,15, 18,57,46); T(i).md(1)=0; % S1.
T(i).Ts(2)=datetime(2018,06,16, 11,55,56); T(i).Te(2)=datetime(2018,06,25, 10,47,12); T(i).md(2)=0; % S2.
T(i).Ts(3)=datetime(2018,06,27, 04,18,53); T(i).Te(3)=datetime(2018,06,30, 08,38,55); T(i).md(3)=0; % S3.
T(i).Ts(4)=datetime(2018,07,01, 04,40,24); T(i).Te(4)=datetime(2018,07,12, 00,50,00); T(i).md(4)=0; % S4.
T(i).Ts(5)=datetime(2018,07,12, 10,24,24); T(i).Te(5)=datetime(2018,07,22, 15,53,06); T(i).md(5)=0; % S5.
T(i).Tf( 1)=datetime(2018,08,21, 01,00,00);
T(i).dtd=0;
T(i).S=[1:5];
i=i+1;

% Helsinki St1 2020.
T(i).Name='St1-2020';
T(i).Filename='/Users/rschultz/Desktop/papers/old/Trailing/IS-bath/data/processed/St1-2020.mat';
T(i).Mc=-1.3;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(2020,05,05, 06,57,44); T(i).Te(1)=datetime(2020,05,21, 03,19,34); T(i).md(1)=0; % S1.
T(i).Tf( 1)=datetime(2020,06,16, 01,00,00);
T(i).dtd=0;
T(i).S=[1];
i=i+1;

% UK-PNR1z.
T(i).Name='PNR1z';
T(i).Filename='/Users/rschultz/Desktop/papers/old/Trailing/IS-bath/data/processed/PNR1z.mat';
T(i).Mc=-0.7;
T(i).Mk=-0.5;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts( 1)=datetime(2018,10,16, 10,26,24); T(i).Te( 1)=datetime(2018,10,16, 14,06,09); T(i).md( 1)=mean([3374 3376]); % S1.
T(i).Ts( 2)=datetime(2018,10,17, 08,12,49); T(i).Te( 2)=datetime(2018,10,17, 16,35,51); T(i).md( 2)=mean([3356 3358]); % S2.
T(i).Ts( 3)=datetime(2018,10,18, 09,55,20); T(i).Te( 3)=datetime(2018,10,18, 14,21,58); T(i).md( 3)=mean([3338 3340]); % S3.
T(i).Ts( 4)=datetime(2018,10,19, 11,32,45); T(i).Te( 4)=datetime(2018,10,19, 14,39,35); T(i).md( 4)=mean([3180 3182]); % S12.
T(i).Ts( 5)=datetime(2018,10,20, 08,06,18); T(i).Te( 5)=datetime(2018,10,20, 10,40,12); T(i).md( 5)=mean([3180 3182]); % S12.
T(i).Ts( 6)=datetime(2018,10,22, 07,46,22); T(i).Te( 6)=datetime(2018,10,22, 11,53,45); T(i).md( 6)=mean([3162 3164]); % S13.
T(i).Ts( 7)=datetime(2018,10,23, 13,32,55); T(i).Te( 7)=datetime(2018,10,23, 14,58,33); T(i).md( 7)=mean([3144 3146]); % S14.
T(i).Ts( 8)=datetime(2018,10,24, 08,09,04); T(i).Te( 8)=datetime(2018,10,24, 12,30,13); T(i).md( 8)=mean([3074 3076]); % S18.
T(i).Ts( 9)=datetime(2018,10,25, 07,39,26); T(i).Te( 9)=datetime(2018,10,25, 11,44,21); T(i).md( 9)=mean([2997 3000]); % S22.
T(i).Ts(10)=datetime(2018,10,25, 14,29,39); T(i).Te(10)=datetime(2018,10,25, 15,42,08); T(i).md(10)=mean([2997 3000]); % S22.
T(i).Ts(11)=datetime(2018,10,26, 07,19,27); T(i).Te(11)=datetime(2018,10,26, 11,41,12); T(i).md(11)=mean([2857 2859]); % S30.
T(i).Ts(12)=datetime(2018,10,27, 08,04,04); T(i).Te(12)=datetime(2018,10,27, 10,53,30); T(i).md(12)=mean([2840 2842]); % S31.
T(i).Ts(13)=datetime(2018,10,29, 08,53,08); T(i).Te(13)=datetime(2018,10,29, 11,39,46); T(i).md(13)=mean([2822 2824]); % S32.
T(i).Ts(14)=datetime(2018,10,30, 08,26,53); T(i).Te(14)=datetime(2018,10,30, 10,07,28); T(i).md(14)=mean([2699 2701]); % S39.
T(i).Ts(15)=datetime(2018,10,31, 09,28,27); T(i).Te(15)=datetime(2018,10,31, 13,34,31); T(i).md(15)=mean([2681 2683]); % S40.
T(i).Ts(16)=datetime(2018,11,02, 08,13,32); T(i).Te(16)=datetime(2018,11,02, 16,43,19); T(i).md(16)=mean([2769 2771]); % S35.
%%%
T(i).Ts(17)=datetime(2018,12,08, 09,24,47); T(i).Te(17)=datetime(2018,12,08, 12,14,03); T(i).md(17)=mean([2734 2736]); % S37.
T(i).Ts(18)=datetime(2018,12,10, 08,34,18); T(i).Te(18)=datetime(2018,12,10, 10,04,58); T(i).md(18)=mean([2734 2736]); % S37.
T(i).Ts(19)=datetime(2018,12,11, 08,32,44); T(i).Te(19)=datetime(2018,12,11, 10,36,32); T(i).md(19)=mean([2716 2718]); % S38.
T(i).Ts(20)=datetime(2018,12,13, 08,30,10); T(i).Te(20)=datetime(2018,12,13, 11,33,42); T(i).md(20)=mean([2699 2701]); % S39.
T(i).Ts(21)=datetime(2018,12,14, 12,27,18); T(i).Te(21)=datetime(2018,12,14, 14,04,44); T(i).md(21)=mean([2681 2683]); % S40.
T(i).Ts(22)=datetime(2018,12,15, 08,59,21); T(i).Te(22)=datetime(2018,12,15, 13,00,13); T(i).md(22)=mean([2663 2665]); % S41.
T(i).Ts(23)=datetime(2018,12,17, 08,13,24); T(i).Te(23)=datetime(2018,12,17, 13,12,47); T(i).md(23)=mean([2663 2665]); % S41.
T(i).Tf( 1)=datetime(2099,01,01, 01,00,00);
T(i).dtd=1;
i=i+1;

% UK-PNR2.
T(i).Name='PNR2';
T(i).Filename='/Users/rschultz/Desktop/papers/old/Trailing/IS-bath/data/processed/PNR2.mat';
T(i).Mc=-1.0;
T(i).Mk=-0.5;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(2019,08,15, 09,39,13); T(i).Te(1)=datetime(2019,08,15, 10,50,48); T(i).md(1)=mean([3106 3111]); % S1.
T(i).Ts(2)=datetime(2019,08,16, 08,01,38); T(i).Te(2)=datetime(2019,08,16, 10,25,14); T(i).md(2)=mean([3092 3097]); % S2.
T(i).Ts(3)=datetime(2019,08,17, 07,59,42); T(i).Te(3)=datetime(2019,08,17, 11,19,49); T(i).md(3)=mean([3077 3082]); % S3.
T(i).Ts(4)=datetime(2019,08,19, 08,12,10); T(i).Te(4)=datetime(2019,08,19, 10,49,39); T(i).md(4)=mean([3063 3068]); % S4.
T(i).Ts(5)=datetime(2019,08,20, 07,57,17); T(i).Te(5)=datetime(2019,08,20, 10,35,01); T(i).md(5)=mean([3048 3053]); % S5.
T(i).Ts(6)=datetime(2019,08,21, 08,36,50); T(i).Te(6)=datetime(2019,08,21, 15,53,28); T(i).md(6)=mean([3033 3038]); % S6.
T(i).Ts(7)=datetime(2019,08,23, 10,21,35); T(i).Te(7)=datetime(2019,08,23, 14,44,56); T(i).md(7)=mean([3019 3024]); % S7.
T(i).Tf( 1)=datetime(2099,01,01, 01,00,00);
T(i).dtd=1;
i=i+1;

% Utah-FORGE 2022.
T(i).Name='FORGE22';
T(i).Filename='data/FORGE2022.mat';
T(i).Mc=-1.2;
T(i).Mk=-1.0;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(2022,04,17, 02,41,36); T(i).Te(1)=datetime(2022,04,17, 05,16,48); T(i).md(1)=mean([10787 10955])*0.3048; % S1.
T(i).Ts(2)=datetime(2022,04,19, 12,50,01); T(i).Te(2)=datetime(2022,04,19, 15,42,50); T(i).md(2)=mean([10560 10580])*0.3048; % S2.
T(i).Ts(3)=datetime(2022,04,21, 13,33,02); T(i).Te(3)=datetime(2022,04,21, 16,19,24); T(i).md(3)=mean([10120 10140])*0.3048; % S3.
T(i).Tf( 1)=datetime(2099,01,01, 01,00,00);
T(i).dtd=1;
i=i+1;

% Utah-FORGE 2024.
T(i).Name='FORGE24';
T(i).Filename='/Users/rschultz/Desktop/papers/old/Trailing/IS-bath/data/processed/FORGE2024.mat';
T(i).Mc=+0.15;
T(i).Mk=+0.30;
T(i).dMb=0.10;
T(i).dMd=0.10;
T(i).Ts( 1)=datetime(2024,04,03, 16,15,06); T(i).Te( 1)=datetime(2024,04,03, 19,56,56); T(i).md( 1)=mean([10120 10140])*0.3048; % S3a(r).
T(i).Ts( 2)=datetime(2024,04,03, 23,56,05); T(i).Te( 2)=datetime(2024,04,04, 03,01,50); T(i).md( 2)=mean([10070 10076])*0.3048; % S4a.
T(i).Ts( 3)=datetime(2024,04,04, 06,47,14); T(i).Te( 3)=datetime(2024,04,04, 09,20,52); T(i).md( 3)=mean([10020 10026])*0.3048; % S5a.
T(i).Ts( 4)=datetime(2024,04,04, 12,31,33); T(i).Te( 4)=datetime(2024,04,04, 22,46,44); T(i).md( 4)=mean([ 9959  9976])*0.3048; % S6a.
T(i).Ts( 5)=datetime(2024,04,05, 06,33,34); T(i).Te( 5)=datetime(2024,04,05, 15,12,12); T(i).md( 5)=mean([ 9798  9901])*0.3048; % S7a.
T(i).Ts( 6)=datetime(2024,04,05, 19,45,07); T(i).Te( 6)=datetime(2024,04,06, 04,01,09); T(i).md( 6)=mean([ 9545  9723])*0.3048; % S8a.
T(i).Ts( 7)=datetime(2024,04,06, 07,00,42); T(i).Te( 7)=datetime(2024,04,06, 13,12,28); T(i).md( 7)=mean([ 9490  9493])*0.3048; % S9a.
T(i).Ts( 8)=datetime(2024,04,07, 02,57,28); T(i).Te( 8)=datetime(2024,04,07, 05,20,44); T(i).md( 8)=mean([ 9270  9323])*0.3048; % S10a.
T(i).Ts( 9)=datetime(2024,04,11, 10,01,52); T(i).Te( 9)=datetime(2024,04,11, 11,12,48); T(i).md( 9)=mean([ ])*0.3048; % S1b.
T(i).Ts(10)=datetime(2024,04,13, 17,10,25); T(i).Te(10)=datetime(2024,04,13, 18,38,43); T(i).md(10)=mean([ ])*0.3048; % S2b.
T(i).Ts(11)=datetime(2024,04,15, 20,24,58); T(i).Te(11)=datetime(2024,04,15, 21,53,12); T(i).md(11)=mean([ ])*0.3048; % S3b.
T(i).Ts(12)=datetime(2024,04,17, 12,02,55); T(i).Te(12)=datetime(2024,04,17, 13,13,02); T(i).md(12)=mean([ ])*0.3048; % S4b.
T(i).Tf( 1)=datetime(2099,01,01, 01,00,00);
T(i).dtd=0;
i=i+1;

% Äspö HRL.
T(i).Name='Aspo';
T(i).Filename='/Users/rschultz/Desktop/papers/old/Trailing/IS-bath/data/processed/Aspo.mat';
T(i).Mc=-4.9;
T(i).dMb=0.1;
T(i).dMd=0.1;
T(i).Ts(1)=datetime(2015,06,03, 09,40,33); T(i).Te(1)=datetime(2015,06,03, 17,00,00); T(i).md(1)=24.80; % HF1.
T(i).Ts(2)=datetime(2015,06,04, 07,15,07); T(i).Te(2)=datetime(2015,06,04, 09,30,00); T(i).md(2)=22.00; % HF2.
T(i).Ts(3)=datetime(2015,06,04, 12,07,02); T(i).Te(3)=datetime(2015,06,04, 14,30,00); T(i).md(3)=19.00; % HF3.
T(i).Ts(4)=datetime(2015,06,09, 10,03,05); T(i).Te(4)=datetime(2015,06,09, 14,00,00); T(i).md(4)=13.65; % HF4.
T(i).Ts(5)=datetime(2015,06,10, 10,35,10); T(i).Te(5)=datetime(2015,06,10, 14,00,00); T(i).md(5)=11.80; % HF5.
T(i).Ts(6)=datetime(2015,06,11, 09,29,03); T(i).Te(6)=datetime(2015,06,12, 12,00,00); T(i).md(6)=04.68; % HF6.
T(i).Tf( 1)=datetime(2015,07,01, 01,01,01);
T(i).dtd=0;
T(i).S=[1];
i=i+1;












