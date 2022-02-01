// Read profile written by vmcw program (x, sx,sy,sz, ntx,nty,ntz, dsx,dsy,dsz,...)
// Write energy parameters: F_SO, F_GR density, mass density
// TODO: move to the main program, use correct units.
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include "grad.h"


int
main(){
  std::cout << "## x  Mt Mn FSO FG1 FG2 FG\n";
  while (std::cin){
    std::string s;
    getline(std::cin, s);
    if (s.size()==0 || s[0]=='#') continue;
    std::istringstream str(s);
    double x, S[3], NT[3], dS[3], dNT[3];
    str >> x >> S[0]   >> S[1]   >> S[2]
             >> NT[0]  >> NT[1]  >> NT[2]
             >> dS[0]  >> dS[1]  >> dS[2]
             >> dNT[0] >> dNT[1] >> dNT[2];
    if (!str) {
      std::cerr << "bad string: " << s << "\n";
      continue;
    }

    double T  = sqrt(pow(NT[0],2) + pow(NT[1],2) + pow(NT[2],2));
    double dT = 0;
    double N[3] = {1,0,0};
    double dN[3] = {0,0,0};
    if (T!=0) {
      for (int i=0;i<3;i++){
        dT += NT[i]*dNT[i]/T;
        N[i] = NT[i]/T;
      }
      for (int i=0;i<3;i++){
        dN[i] = dNT[i]/T - NT[i]*dT/T/T;
      }
    }

    // Mt = chi_B/gamma^2 (th')^2
    // Mn = chi_B/gamma^2 2(1-ct)(n')^2
    double Mt = dT*dT;
    double Mn = 2*(1-cos(T))*(dN[0]*dN[0]+dN[1]*dN[1]+dN[2]*dN[2]);

    // FSO = g_D Delta^2 * 8(ct+0.25)^2
    double FSO = 8*pow(cos(T) + 0.25, 2);

    // FGR = K1/2 Delta^2 FG1 + (K2+K3)/2 Delta^2 FG2 ~= K1 Delta^2 * (FG1/2 + FG2)
    double FG1,FG2;
    fill_eg1_nt_(&FG1, &FG2, N, &T, dN, &dT);

    std::cout << x << " " << Mt << " " << Mn << " " << FSO << " " << FG1 << " " << FG2 << " " << FG1/2.0 + FG2 << "\n";
  }

}