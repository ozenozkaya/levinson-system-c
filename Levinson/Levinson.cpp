#include <iostream>
#include <math.h>

using namespace std;

typedef struct
{
 double daLpcCoeff[11];
 double dGain; //sigma        
}tSLevinsonOutput;

typedef struct
{
        double daAcfCoeff[11];
        unsigned char bPredictor;
        unsigned int  sWindowSample;
}tSLevinsonInput;

tSLevinsonOutput LevinsonFunction(tSLevinsonInput SLevinsonInput)
{
tSLevinsonOutput SLevinsonOutput={0};
unsigned short sCnt=0;
unsigned short sSubCnt=0;
double da_K[10] = {0};
double da_A[10][10]= {0};
double da_E[11] = {0};
double dSum= 0.0;
double dSigmaSquare=0.0;

       da_K[0] = SLevinsonInput.daAcfCoeff[1] / SLevinsonInput.daAcfCoeff[0];   
       da_A[0][0] = da_K[0];
       da_E[1] = (1.0-(da_K[0]*da_K[0]))*SLevinsonInput.daAcfCoeff[0];
       
       for(sCnt=2;sCnt<(unsigned short)SLevinsonInput.bPredictor+1; sCnt++)
       {
          dSum=0;
          for(sSubCnt=1;sSubCnt<(sCnt+1);sSubCnt++)
          {
           dSum=dSum + (da_A[sCnt-2][sSubCnt-1]*SLevinsonInput.daAcfCoeff[sCnt-sSubCnt]);
          }
          da_K[sCnt-1] = (SLevinsonInput.daAcfCoeff[sCnt]-dSum)/da_E[sCnt-1];
          da_A[sCnt-1][sCnt-1] = da_K[sCnt-1];
          for(sSubCnt=1;sSubCnt<(sCnt);sSubCnt++)
          {
           da_A[sCnt-1][sSubCnt-1]=da_A[sCnt-2][sSubCnt-1]-(da_K[sCnt-1]*da_A[sCnt-2][sCnt-sSubCnt-1]);
          }
          da_E[sCnt] = (1.0-(da_K[sCnt-1]*da_K[sCnt-1]))*da_E[sCnt-1];                             
       }
	   SLevinsonOutput.daLpcCoeff[0] = 1.0000;
       for(sSubCnt=0;sSubCnt<(sCnt);sSubCnt++)
       {
        SLevinsonOutput.daLpcCoeff[sSubCnt+1] = (-1.0)*da_A[(unsigned char)SLevinsonInput.bPredictor-1][sSubCnt];
       }
       if(da_E[(unsigned char)SLevinsonInput.bPredictor]>=0.0)
           dSigmaSquare = da_E[(unsigned char)SLevinsonInput.bPredictor];
       else
           dSigmaSquare = (-1.0)*da_E[(unsigned char)SLevinsonInput.bPredictor];
       
       SLevinsonOutput.dGain = sqrt(dSigmaSquare);
       SLevinsonOutput.dGain /= sqrt((float)SLevinsonInput.sWindowSample);    
                 
return SLevinsonOutput;
}

int main()
{
    tSLevinsonInput SMyLevinsonInput={0};
    tSLevinsonOutput SMyLevinsonOutput={0};
    
    SMyLevinsonInput.daAcfCoeff[0]=0.0011784;
    SMyLevinsonInput.daAcfCoeff[1]=0.0011253;
    SMyLevinsonInput.daAcfCoeff[2]=0.0010583;
    SMyLevinsonInput.daAcfCoeff[3]=0.0009857;
    SMyLevinsonInput.daAcfCoeff[4]=0.0009287;
    SMyLevinsonInput.daAcfCoeff[5]=0.0008846;
    SMyLevinsonInput.daAcfCoeff[6]=0.0008563;
    SMyLevinsonInput.daAcfCoeff[7]=0.0008332;
    SMyLevinsonInput.daAcfCoeff[8]=0.0008079;
    SMyLevinsonInput.daAcfCoeff[9]=0.0007686;
    SMyLevinsonInput.daAcfCoeff[10]=0.0007303;
    
    SMyLevinsonInput.bPredictor = 10;
    SMyLevinsonInput.sWindowSample = 240;
    
    SMyLevinsonOutput = LevinsonFunction(SMyLevinsonInput);
    printf("Kazanc: %f\n\r",SMyLevinsonOutput.dGain);
    for(int sCnt=0;sCnt<11;sCnt++)
           printf("%d. Katsayi: %f\n\r",sCnt,SMyLevinsonOutput.daLpcCoeff[sCnt]); 
    system("PAUSE");
    return 0;
    }
