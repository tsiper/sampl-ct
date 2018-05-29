The following .m files are given as examples:

SPURS_Spiral20k.m - Performs reconstruction from 20000 samples of an analyitcal brain phantom taken on an Arcemidain spiral trajectory with ISNR = 30dB.
SPURS_Spiral30k.m - Performs reconstruction from 30000 samples of an analyitcal brain phantom taken on an Arcemidain spiral trajectory with ISNR = 30dB.
SPURS_Spiral50k.m - Performs reconstruction from 50000 samples of an analyitcal brain phantom taken on an Arcemidain spiral trajectory with ISNR = 30dB.
SPURS_Radial30k.m - Performs reconstruction from 30661 samples of an analyitcal brain phantom taken on a radial trajectory with 60 spokes and with ISNR = 30dB.

The input SNR can be changed from 30 to any value desired (noiseless = inf).
The input pahntom can be chaned to from 'Brain' to a box phantom - 'Phantom1' or to an analytical Shep loagen pahntom - 'AnalyticalSL'. 
The degree of the spline funciton is controlled by the values of  - BsplineDegree.
The value fo the regularization parameter is - Rho.
The number of iteration is - Niterations.
The over grid factor is - OverGridFactor.

Notes:
* In the first run the PHI matrix for each trajectory is calculated. This is a time consuming task.
* In order to save time for samples taken on the same trajectory the LU factor and PHI matrix are saved on disk (in \PreCalculated\LU and \PreCalculated\PHI folders). 
* Once the simulation re-run with the same trajectory these matrices are reloaded from file and reused. You can disable this option or delete the generated files if you want and the files will be calculated again when rquired.
* The same is done for phantom data for a given trajectory (\PreCalculated\PhantomData). 
The 'Brain' and 'AnalyticalSL' phantom data is generated using http://bigwww.epfl.ch/algorithms/mriphantom/ from EPFL.
This demo is provided with samples for the four trajectories (as .mat files). 
If a different trajectory is required you can either use the 'Phantom1' phantom or download the package (make sure to generate data for a single coil with uniform sensitivity).
* Matlab 2014a or newer is required.
