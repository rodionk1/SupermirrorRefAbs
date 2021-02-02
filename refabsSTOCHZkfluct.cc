#include <fstream>
#include <iostream>
#include <math.h>
#include <iomanip>
#include <complex>
#include <cmath>
#include <string>
#include <random>
#include <omp.h>
#include <chrono>

#define MMIN 0.8
#define MMAX 1.5
#define NSIM 500 // Number of simulation for layer thickness randomization
#define NTKSIM 50 // iterations for k and thickness fluctuations
#define DEPACCUR 0.005L // accuracy in layer thickness deposition
#define ISIGMA 2.0L // Roughness of the substrate, in AA
#define ISIGMAPLUS 1.0E-5 // Roughness increase in AA per AA (1.0e-4 = 1 AA/mkm)
#define KSPREADREL 0.00L // fluctuations in kperp relative to mean value
#define KSPREADABS 0.03L*0.0109L // fluctuations in kperp, absolute, in AA-1
#define KSPREADABS 0.00L // fluctuations in kperp, absolute, in AA-1
#define LAMBDA 5.0L // Wavelength in AA
#define LAMBDA2200 1.798L // Wavelength of thermal neutrons in AA
#define SIGSNI 18.5 //Elastic scattering 5AA
#define SIGSMO 5.71
#define SIGSTI 4.35

#define SIGDNIMO 0.453 // cm-1, macroscopic diffuse scattering cross section in NiMo layer
#define SIGDTI 0.163 // same for Ti
#define SIGANINIMO 0.382 // cm-1, macroscopic absorption cross section for Ni in NiMo layer
#define SIGAMONIMO 0.023 // cm-1, macroscopic absorption cross section for Mo in NiMo layer
#define SIGATI 0.346 // same for Ti
/*Elastic scattering X-sections: ENDF
L    Ni   Mo   Ti   
5AA 19.07 5.21 4.73
4AA 18.63 5.14 4.60
3AA 18.31 5.08 4.51
2AA 18.06 5.03 4.46
1AA 17.93 5.01 4.41
Sears:
Thermal 18.5 5.71 4.35
*/

#define PI 3.14159265L // Pi value, recognized as long double

std::complex<double> qform (std::complex<double> v1[2], std::complex<double> A[2][2], std::complex<double> v2[2])
{
  std::complex<double> result (0.0,0.0);
  for (int i=0;i<2;i++){
    for (int j=0; j<2;j++){
      result+= v1[i]*A[i][j]*v2[j];
    }
  }
  return result;
}




void matmatmult (std::complex<long double> A[2][2], std::complex<long double> B[2][2], std::complex<long double> R[2][2])
{
  std::complex<long double> li(0.0,1.0);
  for (int i=0;i<2;i++){
    for (int j=0; j<2;j++){
      long double zero=0.0;
     R[i][j]= zero*li;
      for (int m=0; m<2; m++){
	R[i][j]+=A[i][m]*B[m][j];
      }
    }
  }
}



void matvecmult (std::complex<long double> A[2][2], std::complex<long double> V[2], std::complex<long double> VR[2])
{
  std::complex<long double> li(0.0,1.0);
  for (int i=0;i<2;i++){
       VR[i]= 0.0L+0.0L*li;
      for (int m=0; m<2; m++){
	VR[i]+=A[i][m]*V[m];
      }
  }
}



int main (int argc, char* argv[])
{
  std::random_device rd;
  typedef std::chrono::high_resolution_clock myclock;
  myclock::time_point beginning = myclock::now();
  //  std::mt19937 gen(rd());
  std::normal_distribution<long double> gauss(0.0L,1.0L);
  std::complex<long double> il(0.0,1.0);
  //Reading layer parameters from file 
  std::ifstream layers;
  layers.open(argv[1]);
  //First line -- number of layers and substrate qc.

  int Nlayers;
  double qc_sub;
  
  layers >> Nlayers >> qc_sub;
  
  long double* thicknessfile = new long double [Nlayers+2];
  long double* roughness = new long double [Nlayers+2];
  long double* qc2 = new long double [Nlayers+2];
  long double* Sigma_s = new long double [Nlayers+2];
  long double* Sigma_a = new long double [Nlayers+2];
  long double* rN = new long double [Nlayers+2];
  
  long double* CurrentAV = new long double [Nlayers+2];
  long double* CurrentOUTAV = new long double [Nlayers+2];
  

      long double** Current  = new long double*[NTKSIM];
      for (int i=0; i < NTKSIM; i++){
	Current[i]=new long double[Nlayers+2];}

      long double** CurrentOUT  = new long double*[NTKSIM];
      for (int i=0; i < NTKSIM; i++){
	CurrentOUT[i]=new long double[Nlayers+2];}

  
  int Num_layer; // layer number
  long double d_layer;   // layer thickness
  std::string s1,s2; // just string
  char type_layer; 
    long double Sig_s_M=SIGDNIMO*1.0E-08;//1.6346E-08; // in AA-1
    long double Sig_s_V=SIGDTI*1.0E-08;// in AA-1
 
  roughness[Nlayers+1] = ISIGMA; // roughness of the substrate
  while (layers >> Num_layer >> d_layer >> s1 >> type_layer >> s2)
    {
      thicknessfile[Nlayers+1-Num_layer] = d_layer;
      roughness[Nlayers+1-Num_layer] = roughness[Nlayers+2-Num_layer] + d_layer*ISIGMAPLUS;
      //    std::cout << Num_layer <<" " << d_layer <<  " " << type_layer << std::endl;
      //      std::cout << Nlayers+1-Num_layer << " " <<  roughness[Nlayers+1-Num_layer] << std::endl;
      if (type_layer=='M')
	{ // if NiMo layer
	  qc2[Nlayers+1-Num_layer] = 4.0*PI* 0.93951E-05; // in AA-2
	  Sigma_s[Nlayers+1-Num_layer] = Sig_s_M;
	  // Sigma_a[Nlayers+1-Num_layer] = 4.0582E-09*LAMBDA/LAMBDA2200;
	  Sigma_a[Nlayers+1-Num_layer] = (SIGANINIMO+SIGAMONIMO)*1.E-08*LAMBDA/LAMBDA2200;
	  rN[Nlayers+1-Num_layer] =   0.93951E-05; // in AA-2       
	};
      if (type_layer=='V')
	{ // if Ti layer
	  qc2[Nlayers+1-Num_layer] = -4.0*PI* 0.19524E-05; // in AA-2
	  Sigma_s[Nlayers+1-Num_layer] = Sig_s_V;//0.24733E-08; // in AA-1
	  Sigma_a[Nlayers+1-Num_layer] = SIGATI*1.E-08*LAMBDA/LAMBDA2200;
	  rN[Nlayers+1-Num_layer] =   -0.19524E-05; // in AA-2       
	};
    };
  roughness[0]=0.0;
  qc2[0]=0.0; // vacuum
  thicknessfile[0]=0.0;
  Sigma_s[0] = 0.0; // vacuum
  Sigma_a[0]=0.0; // vacuum
  rN[0] = 0.0; // vacuum
  
  qc2[Nlayers+1] = qc_sub*qc_sub; // shift in the substrate

  thicknessfile[Nlayers+1]=0.0;
  Sigma_s[Nlayers+1] = 0.0; // will not be used, neglect in the first approximation
  Sigma_a[Nlayers+1]=4.0E-08; // absorption cross section for boron in Borkron glass is approx 3.9cm-1 = 3.9E-08 AA-1
  rN[Nlayers+1] = 0.0; // will not be used in the first approximation


  long double kperp; // transverse part of wave vector in vacuum

std::cout << "q | m | R | Loss_NiMo | %Loss_NiMo | Loss_Ti | %Loss_Ti | Loss_transmitted | %transmitted | Scatt_roughnes | %capture in Ni | %capture in Mo | %capture in Ti | abs_Ni per incident | abs_Mo per incident | abs_Ti per incident" <<std::endl;


  for (int iq=200*MMIN; iq<=200*MMAX; iq++) // cycle in transverse wave vector
    {
      kperp = 0.00109+0.05*0.00109*iq; // transverse wave vector in AA-1

      for (int i=0; i<=Nlayers+1; i++){
      	for (int j =0; j<NTKSIM; j++){
	  Current[j][i] = 0.0;
	  CurrentOUT[j][i]=0.0;
	};
      };

      //setting values of the matrix      
      long double value [NTKSIM];
      long double LossNiMo = 0.0L;
      long double LossTi = 0.0L;

	myclock::duration d = myclock::now() - beginning;
	unsigned seed = d.count();
         std::mt19937 gen(seed); 
 

      long double ksim[NTKSIM];
      long double** tsim  = new long double*[NTKSIM];
      for (int i=0; i < NTKSIM; i++){
	tsim[i]=new long double[Nlayers+2];}


      for (int iksim=0; iksim < NTKSIM; iksim++){
      ksim[iksim]=kperp+KSPREADREL*kperp*gauss(gen)+KSPREADABS*gauss(gen); // introducing spread in transverse momentum 
  for (int i=Nlayers; i>0; i--)  tsim[iksim][i]= thicknessfile[i]+DEPACCUR*gauss(gen)*thicknessfile[i];
}

#pragma omp parallel shared(Current,CurrentOUT,CurrentOUTAV, CurrentAV, value) 
      {
#pragma omp for schedule(dynamic,1) nowait
  for (int iksim=0; iksim < NTKSIM; iksim++){
   std::complex<long double> T[2][2]; // T-matrix for multilayer
      std::complex<long double> temp[2][2]; // for storing intermediate results of the multiplication
      std::complex<long double> TR[2][2]; // T-matrix for multilayer
      std::complex<long double> tempR[2][2]; // for storing intermediate results of the multiplication


  std::complex <long double>* Vup = new std::complex <long double> [Nlayers+2];
  std::complex <long double>* Vdown = new std::complex <long double> [Nlayers+2];
  std::complex <long double>* Uup = new std::complex <long double> [Nlayers+2];
  std::complex <long double>* Udown = new std::complex <long double> [Nlayers+2];
  std::complex <long double>* alpha = new std::complex <long double> [Nlayers+2];
  std::complex <long double>* beta = new std::complex <long double> [Nlayers+2];


  long double* thickness = new long double [Nlayers+2];
  long double* zinterface = new long double [Nlayers+2];
  long double* zinterfacesim = new long double [Nlayers+2];
  thickness[Nlayers+1] = 0.0; // 1 cm, but set tozero to get right substrate position variation in the calculation
  thickness[0]=0.0; // set to zero by convention, see Shelten and Mika
  zinterface[Nlayers+1]=0; // start counting coating thickness from the substrate surface
      for (int i=0; i<=Nlayers+1; i++){
	alpha[i]=0.0;
	beta[i]= 0.0;
      };

 // qperp in the layers
  std::complex <long double> * qperp = new std::complex <long double> [Nlayers+2];
  //T-matrices for the layers  
  std::complex <long double>  Tj[2][2];
  //  std::complex <long double>  TjR[2][2];


      // Wave vectors in the layers (complex values)

  qperp[0]= ksim[iksim];
      for (int i=1; i<Nlayers+2; i++){
	qperp[i] =sqrt(qperp[0]*qperp[0] - qc2[i]*1.0L + il*2.0L*PI/LAMBDA*(Sigma_s[i]+Sigma_a[i])*(1.0L-0.5L*LAMBDA*LAMBDA*rN[i]/3.14159265L));
     
      };

      for (int i=0; i<=Nlayers+1; i++){
	alpha[i]=0.0;
	beta[i]= 0.0;
      };      
      alpha[0]=1.0;

 for (int i=Nlayers+1; i>0; i--){
    zinterface[i] = 0.0;
  }
 

  for (int i=Nlayers; i>0; i--){
    thickness[i] = tsim[iksim][i];
    zinterface[i] = zinterface[i+1]+thickness[i];
  }
  zinterface[0]= zinterface[1];
    for (int isim=0; isim< NSIM; isim++){

     
      T[0][0]=1;
      T[1][1]=1;
      T[1][0]=0;
      T[0][1]=0;

      TR[0][0]=1;
      TR[1][1]=1;
      TR[1][0]=0;
      TR[0][1]=0;


      long double onehalf = 0.5;
      Vup[0]=1.0L;
      Vdown[0]=0.0L;
      Uup[0]=0.0L;
      Udown[0]=1.0L;
	myclock::duration d1 = myclock::now() - beginning;
	unsigned seed1 = d1.count();
         std::mt19937 gen1(seed1); 


      //generating a sequence of interface boundaries
        for (int i = Nlayers+1; i>0; i--){
	  zinterfacesim[i] = zinterface[i]+roughness[i]*gauss(gen1);
	}
	//       zinterfacesim[0]=zinterfacesim[1];
   zinterfacesim[0]=zinterface[1]; // any fixed number in fact, just to get the plane wave shift when it touches the rough surface
       // double totalthick = 0.0;
      for (int i=1; i<Nlayers+2; i++){
	long double lthickness = zinterfacesim[i-1] - zinterfacesim[i];
	Tj[0][0] = (qperp[i]+qperp[i-1])/qperp[i]*onehalf*exp(il*qperp[i-1]*lthickness);
	Tj[0][1] = (qperp[i]-qperp[i-1])/qperp[i]*onehalf*exp(-il*qperp[i-1]*lthickness);
	Tj[1][0] = (qperp[i]-qperp[i-1])/qperp[i]*onehalf*exp(il*qperp[i-1]*lthickness);
	Tj[1][1] = (qperp[i]+qperp[i-1])/qperp[i]*onehalf*exp(-il*qperp[i-1]*lthickness);
	
// compensating precision loss by interrupting layer sequence at few first layers for small momentum transfer
	if((real(qperp[0])<0.00109L*6.5L)&&(i>7)){ 
	  Vup[i]=0.0;
	  Vdown[i]=0.0;
	  Uup[i]=0.0;
	  Udown[i]=0.0;
	}else {
	temp[0][0]=T[0][0];
	temp[1][1]=T[1][1];
	temp[1][0]=T[1][0];
	temp[0][1]=T[0][1];
	matmatmult(Tj,temp,T);

	std::complex <long double> Vstart[2]={1.0L,0.0L};
	std::complex <long double> Ustart[2]={0.0L,1.0L};

	std::complex <long double> Vend[2];
	std::complex <long double> Uend[2];

	matvecmult(T,Vstart,Vend);
	matvecmult(T,Ustart,Uend);
	
	Vup[i]=Vend[0];
	Vdown[i]=Vend[1];
	Uup[i]=Uend[0];
	Udown[i]=Uend[1];
	}
      }


      
      std::complex <long double> beta0 = -T[1][0]/T[1][1];

      // alpha0 =1
      //Now calculating coefficients alpha and beta
      beta[0]+=beta0;
      for (int i=1; i<=Nlayers+1; i++){
	alpha[i]+=(1.0L*Vup[i]+beta0*Uup[i])*exp(-il*qperp[i]*(zinterfacesim[i]-zinterface[i]));
	beta[i]+=(1.0L*Vdown[i]+beta0*Udown[i])*exp(il*qperp[i]*(zinterfacesim[i]-zinterface[i]));
      };
   }//isim
 

      //Now calculating currents at the layer interfaces
    beta[0]/=NSIM;
      for (int i=1; i<=Nlayers+1; i++){
	alpha[i]/=NSIM;
	beta[i]/=NSIM;
      	Current[iksim][i] = real( qperp[i]*(alpha[i]-beta[i])*(conj(alpha[i])+conj(beta[i])))/real(qperp[0])/NTKSIM;
      	CurrentOUT[iksim][i] = real( qperp[i]*(alpha[i]*exp(il*qperp[i]*thickness[i])-beta[i]*exp(-il*qperp[i]*thickness[i]))*(conj(alpha[i]*exp(il*qperp[i]*thickness[i]))+conj(beta[i]*exp(-il*qperp[i]*thickness[i]))))/real(qperp[0])/NTKSIM;
      };

  value[iksim] = norm(beta[0])/NTKSIM; // this is reflectivity

  delete[] qperp;
  delete[] Vup;
  delete[] Vdown;
  delete[] Uup;
  delete[] Udown;
  delete[] alpha;
  delete[] beta;
  delete[] thickness;
  delete[] zinterface;
  delete[] zinterfacesim;
      }//iksim
   } // pragma omp parallel

      long double Refl = 0.0;
      for (int iksim=0; iksim< NTKSIM; iksim++) Refl+=value[iksim];
      
      for(int j=1; j<=Nlayers+1; j++){
	CurrentOUTAV[j]=0;
	CurrentAV[j]=0;
	for (int i=0; i<NTKSIM; i++) {
	  CurrentOUTAV[j]+=CurrentOUT[i][j];
	  CurrentAV[j]+=Current[i][j];	
	    }
      }
      for (int i=0; i<NTKSIM; i++) delete[] tsim[i];
      
     delete[] tsim;
      
      for (int i=1; i<=(Nlayers+1)/2; i++){
	LossNiMo+=CurrentAV[2*i-1]-CurrentOUTAV[2*i-1];
      }

      for (int i=1; i<=Nlayers/2; i++){
	LossTi+=CurrentAV[2*i]-CurrentOUTAV[2*i];
      }

      long double absfracNi = SIGANINIMO*LAMBDA/LAMBDA2200/(SIGDNIMO + (SIGANINIMO+SIGAMONIMO)*LAMBDA/LAMBDA2200);
      long double absfracMo = SIGAMONIMO*LAMBDA/LAMBDA2200/(SIGDNIMO + (SIGANINIMO+SIGAMONIMO)*LAMBDA/LAMBDA2200);
      long double absfracTi = SIGATI*LAMBDA/LAMBDA2200/(SIGDTI + SIGATI*LAMBDA/LAMBDA2200);
      //      std::cout << "q | m | R | Loss_NiMo | %Loss_NiMo | Loss_Ti | %Loss_Ti | Loss_transmitted | %transmitted | Scatt_roughnes | %capture in Ni | %capture in Mo | %capture in Ti | abs_Ni per incident | abs_Mo per incident | abs_Ti per incident" <<std::endl;

      

std::cout << 2.0*kperp<<" " << kperp/0.0109<< " " << Refl <<" " 
	  << LossNiMo <<" " << LossNiMo/(1-Refl) <<" " << LossTi <<" " << LossTi/(1-Refl) << " "
		<< CurrentAV[Nlayers+1] <<  " " << CurrentAV[Nlayers+1]/(1-Refl) <<" " << 1.0-(Refl+LossNiMo+LossTi+CurrentAV[Nlayers+1]) 
		<< " " << absfracNi*LossNiMo/(1-Refl) << " " << absfracMo*LossNiMo/(1-Refl) << " " << absfracTi*LossTi/(1-Refl)
		<< " " << absfracNi*LossNiMo << " " << absfracMo*LossNiMo<< " " << absfracTi*LossTi;
	std::cout<< std::endl;
    } // end of iteration in momentum
  layers.close();

      for (int i=0; i<NTKSIM; i++) delete[] Current[i];      
     delete[] Current;

      for (int i=0; i<NTKSIM; i++) delete[] CurrentOUT[i];      
     delete[] CurrentOUT;

  delete[] roughness;
  delete[] qc2;
  delete[] Sigma_s;
  delete[] Sigma_a;
  delete[] rN;

  
  delete[] CurrentAV;
  delete[] CurrentOUTAV;
}
