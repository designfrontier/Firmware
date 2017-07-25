/*This file is part of the Maslow Control Software.

    The Maslow Control Software is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Maslow Control Software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Maslow Control Software.  If not, see <http://www.gnu.org/licenses/>.

    Copyright 2014-2017 Bar Smith*/

/*
The Kinematics module relates the lengths of the chains to the position of the cutting head
in X-Y space.
*/

#include "Arduino.h"
#include "Kinematics.h"


Kinematics::Kinematics(){

}

void Kinematics::_verifyValidTarget(float* xTarget,float* yTarget){
    //If the target point is beyond one of the edges of the board, the machine stops at the edge

    *xTarget = (*xTarget < -halfWidth) ? -halfWidth : (*xTarget > halfWidth) ? halfWidth : *xTarget;
    *yTarget = (*yTarget < -halfHeight) ? -halfHeight : (*yTarget > halfHeight) ? halfHeight : *yTarget;

}

void Kinematics::recomputeGeometry(){
    /*
    Some variables are computed on class creation for the geometry of the machine to reduce overhead,
    calling this function regenerates those values.
    */
    h = sqrt((l/2)*(l/2) + s * s);
    Theta = atan(2*s/l);
    Psi1 = Theta - Phi;
    Psi2 = Theta + Phi;

}

void  Kinematics::inverse(float xTarget,float yTarget, float* aChainLength, float* bChainLength){

    //Confirm that the coordinates are on the wood
    _verifyValidTarget(&xTarget, &yTarget);
    
    float Chain1 = sqrt(pow((-1*_xCordOfMotor - xTarget),2)+pow((_yCordOfMotor - yTarget),2));
    float Chain2 = sqrt(pow((_xCordOfMotor - xTarget),2)+pow((_yCordOfMotor - yTarget),2));
    
    //Subtract of the virtual length which is added to the chain by the rotation mechanism
    Chain1 = Chain1 - rotationDiskRadius;
    Chain2 = Chain2 - rotationDiskRadius;
    
    *aChainLength = Chain1;
    *bChainLength = Chain2;

}

void  Kinematics::forward(const float& chainALength, const float& chainBLength, float* xPos, float* yPos){

    float xGuess = 0;
    float yGuess = 0;

    float guessLengthA;
    float guessLengthB;

    int guessCount = 0;
    int maxNumberOfGuesses = 200;

    while(1){


        //check our guess
        inverse(xGuess, yGuess, &guessLengthA, &guessLengthB);

        float aChainError = chainALength - guessLengthA;
        float bChainError = chainBLength - guessLengthB;


        //adjust the guess based on the result
        xGuess = xGuess + .1*aChainError - .1*bChainError;
        yGuess = yGuess - .1*aChainError - .1*bChainError;

        guessCount++;

        //Prevent the connection from timing out
        Serial.print(F("[PEk:"));
        Serial.print(aChainError);
        Serial.print(',');
        Serial.print(bChainError);
        Serial.print(',');
        Serial.print('0');
        Serial.println(F("]"));

        //if we've converged on the point...or it's time to give up, exit the loop
        if((abs(aChainError) < .1 && abs(bChainError) < .1) or guessCount > maxNumberOfGuesses){
            if(guessCount > maxNumberOfGuesses){
                Serial.println(F("Message: Unable to find valid machine position. Please calibrate chain lengths."));
                Serial.println(F("Lengths: "));
                Serial.println(chainALength);
                Serial.println(chainBLength);
                *xPos = 0;
                *yPos = 0;
            }
            else{
                *xPos = xGuess;
                *yPos = yGuess;
            }
            break;
        }
    }
}

void  Kinematics::_MatSolv(){
    float Sum;
    int NN;
    int i;
    int ii;
    int J;
    int JJ;
    int K;
    int KK;
    int L;
    int M;
    int N;

    float fact;

    // gaus elimination, no pivot

    N = 3;
    NN = N-1;
    for (i=1;i<=NN;i++){
        J = (N+1-i);
        JJ = (J-1) * N-1;
        L = J-1;
        KK = -1;
        for (K=0;K<L;K++){
            fact = Jac[KK+J]/Jac[JJ+J];
            for (M=1;M<=J;M++){
                Jac[KK + M]= Jac[KK + M] -fact * Jac[JJ+M];
            }
        KK = KK + N;
        Crit[K] = Crit[K] - fact * Crit[J-1];
        }
    }

//Lower triangular matrix solver

    Solution[0] =  Crit[0]/Jac[0];
    ii = N-1;
    for (i=2;i<=N;i++){
        M = i -1;
        Sum = Crit[i-1];
        for (J=1;J<=M;J++){
            Sum = Sum-Jac[ii+J]*Solution[J-1];
        }
    Solution[i-1] = Sum/Jac[ii+i];
    ii = ii + N;
    }
}

float Kinematics::_moment(const float& Y1Plus, const float& Y2Plus, const float& Phi, const float& MSinPhi, const float& MSinPsi1, const float& MCosPsi1, const float& MSinPsi2, const float& MCosPsi2){   //computes net moment about center of mass
    float Temp;
    float Offsetx1;
    float Offsetx2;
    float Offsety1;
    float Offsety2;
    float Psi1;
    float Psi2;
    float TanGamma;
    float TanLambda;

    Psi1 = Theta - Phi;
    Psi2 = Theta + Phi;

    Offsetx1 = h * MCosPsi1;
    Offsetx2 = h * MCosPsi2;
    Offsety1 = h * MSinPsi1;
    Offsety2 = h * MSinPsi2;
    TanGamma = (y - Offsety1 + Y1Plus)/(x - Offsetx1);
    TanLambda = (y - Offsety2 + Y2Plus)/(D -(x + Offsetx2));

    return h3*MSinPhi + (h/(TanLambda+TanGamma))*(MSinPsi2 - MSinPsi1 + (TanGamma*MCosPsi1 - TanLambda * MCosPsi2));
}

void Kinematics::_MyTrig(){
    float Phisq = Phi * Phi;
    float Phicu = Phi * Phisq;
    float Phidel = Phi + DeltaPhi;
    float Phidelsq = Phidel * Phidel;
    float Phidelcu = Phidel * Phidelsq;
    float Psi1sq = Psi1 * Psi1;
    float Psi1cu = Psi1sq * Psi1;
    float Psi2sq = Psi2 * Psi2;
    float Psi2cu = Psi2 * Psi2sq;
    float Psi1del = Psi1 - DeltaPhi;
    float Psi1delsq = Psi1del * Psi1del;
    float Psi1delcu = Psi1del * Psi1delsq;
    float Psi2del = Psi2 + DeltaPhi;
    float Psi2delsq = Psi2del * Psi2del;
    float Psi2delcu = Psi2del * Psi2delsq;

    // Phirange is 0 to -27 degrees
    // sin -0.1616   -0.0021    1.0002   -0.0000 (error < 6e-6)
    // cos(phi): 0.0388   -0.5117    0.0012    1.0000 (error < 3e-5)
    // Psi1 range is 42 to  69 degrees,
    // sin(Psi1):  -0.0942   -0.1368    1.0965   -0.0241 (error < 2.5 e-5)
    // cos(Psi1):  0.1369   -0.6799    0.1077    0.9756  (error < 1.75e-5)
    // Psi2 range is 15 to 42 degrees
    // sin(Psi2): -0.1460   -0.0197    1.0068   -0.0008 (error < 1.5e-5)
    // cos(Psi2):  0.0792   -0.5559    0.0171    0.9981 (error < 2.5e-5)

    MySinPhi = -0.1616*Phicu - 0.0021*Phisq + 1.0002*Phi;
    MySinPhiDelta = -0.1616*Phidelcu - 0.0021*Phidelsq + 1.0002*Phidel;

    SinPsi1 = -0.0942*Psi1cu - 0.1368*Psi1sq + 1.0965*Psi1 - 0.0241;//sinPsi1
    CosPsi1 = 0.1369*Psi1cu - 0.6799*Psi1sq + 0.1077*Psi1 + 0.9756;//cosPsi1
    SinPsi2 = -0.1460*Psi2cu - 0.0197*Psi2sq + 1.0068*Psi2 - 0.0008;//sinPsi2
    CosPsi2 = 0.0792*Psi2cu - 0.5559*Psi2sq + 0.0171*Psi2 + 0.9981;//cosPsi2

    SinPsi1D = -0.0942*Psi1delcu - 0.1368*Psi1delsq + 1.0965*Psi1del - 0.0241;//sinPsi1
    CosPsi1D = 0.1369*Psi1delcu - 0.6799*Psi1delsq + 0.1077*Psi1del + 0.9756;//cosPsi1
    SinPsi2D = -0.1460*Psi2delcu - 0.0197*Psi2delsq + 1.0068*Psi2del - 0.0008;//sinPsi2
    CosPsi2D = 0.0792*Psi2delcu - 0.5559*Psi2delsq + 0.0171*Psi2del +0.9981;//cosPsi2

}

float Kinematics::_YOffsetEqn(const float& YPlus, const float& Denominator, const float& Psi){
    float Temp;
    Temp = ((sqrt(YPlus * YPlus - R * R)/R) - (y + YPlus - h * sin(Psi))/Denominator);
    return Temp;
}
