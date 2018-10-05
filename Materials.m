function Materials;

MantleRhoS = 3400;
MantleRhoM = 3400;
MantleCp   = 1131;
MantleL    = 0.0;
Mantlek    = 3.4;
MantleXTslope = 1.0e-3.*[1 1];
MantleXTinter = -1e3.*[1 1];
MantleMeltCompaction = [0.00];

BasaltRhoS = 3100;
BasaltRhoM = 2830;
BasaltCp   = 1480;
BasaltL    = 4.0e5;
Basaltk    = 2.6;
BasaltXTslope = [0.00178571 0.00454545];
BasaltXTinter = [-1.7913393 -5.4682168];
BasaltMeltCompaction = [0.10];

LCrustRhoS = 3050;
LCrustRhoM = 2300;
LCrustCp   = 1390;
LCrustL    = 3.5e5;
LCrustk    = 2.6;
LCrustXTslope = 0.0012.*[1 1];
LCrustXTinter = -1.2300.*[1 1];
LCMeltCompaction = [0.00];

UCrustRhoS = 2650;
UCrustRhoM = 2650;
UCrustCp   = 1370;
UCrustL    = 2.7e5;
UCrustk    = 3.0;
UCrustXTslope = [0.004 0.002631579];
UCrustXTinter = [-3.8126 -2.337236842];
UCMeltCompaction = [0.00];

save MaterialLibrary.mat

clear
