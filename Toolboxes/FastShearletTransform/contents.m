% FFST m-files
%
%--------------------------------------------------------------------------
% Sören Häuser ~ FFST ~ 2014-07-22 ~ last edited: 2014-07-22 (Sören Häuser)

% FFST/
% |- create/
% |  |- myBall.m                          - create a ball/circle
% |  |- myPicture.m                       - create an image with different shapes
% |  |- myPicture2.m                      - create an image with other shapes
% |  |- myRhombus.m                       - create a rhombus
% |  |- mySquare.m                        - create a square
% |- helper/
% |  |- checkInputs                       - validate inputs
% |  |  |- checkCoefficients.m            - 
% |  |  |- checkImage.m                   - 
% |  |  |- checkLength.m                  - 
% |  |  |- checkNumOfScales.m             -
% |  |  |- checkShearletSpect.m           -
% |  |  |- defaultNumberOfScales.m        - compute number of scales
% |  |- parseShearletParameterInputs.m    - parse Inputs
% |  |- scalesShearsAndSpectra.m          - compute the spectra of the shearlets for all a and s
% |  |- shearletScaleShear.m              - convert index to a and s or j and k and vice versa
% |- shearlets/
% |  |- bump.m                            - psi_2
% |  |- meyeraux.m                        - v
% |  |- meyerNewShearletSpect.m           - smooth shearlet
% |  |- meyerScaling.m                    - phi
% |  |- meyerShearletSpect.m              - default shearlet
% |  |- meyerWavelet.m                    - psi_1
% | 
% |- inverseShearletTransformSpect.m      - inverse shearlet transform
% |- shearletTransformSpect.m             - shearlet transform


