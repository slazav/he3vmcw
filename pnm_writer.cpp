#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "pnm_writer.h"

// convert 3D unit vector to a color
int vec_to_cal(const double vx, const double vy, const double vz){
  double a = 3*(atan2(vx,vy)/M_PI+1.0); // azimuth: 0..6
  int ai = floor(a); // color range 0..5
  int d = floor((a-ai)*256);
  int r=255,g=0,b=0;
  switch (ai) {
    case 0: r=255; g = d; b=0; break;
    case 1: r=255-d; g=255; b=0; break;
    case 2: r=0; g=255; b=floor(d); break;
    case 3: r=0; g=255-d; b=255; break;
    case 4: r=d; g=0; b=255; break;
    case 5: r=255; g=0; b=255-d; break;
  }
  double f = 2*atan2(vz, hypot(vx,vy))/M_PI; // -1..1
  if (f<0) {
    r = r*(1.0+f);
    g = g*(1.0+f);
    b = b*(1.0+f);
  }
  if (f>0) {
    r = 255.0*f + r*(1.0-f);
    g = 255.0*f + g*(1.0-f);
    b = 255.0*f + b*(1.0-f);
  }

  if (r<0) r=0; if (r>255) r=255;
  if (g<0) g=0; if (g>255) g=255;
  if (b<0) b=0; if (b>255) b=255;
  return (r<<16) + (g<<8) + b;
}

pnm_writer::pnm_writer(const char *fname_):H(0),W(0),fname(fname_){
}

void
pnm_writer::write(const std::vector<double> zsol,
                  const std::vector<double> usol, int NPDE){

  // open file for appending; write to its end
  std::ofstream ss;
  std::ios_base::openmode flags =
     (H==0)? std::fstream::out | std::fstream::trunc :
             std::fstream::out | std::fstream::in | std::fstream::ate;
  ss.open (fname, flags);

  // if it is the first line
  if (H==0){
    ss.seekp(0);
    W = zsol.size()*2+1;
    ss << "P6\n" << W << " ";
    width_pos = ss.tellp();
    ss << "                  \n255\n";
  }

  int x0=100, y0=100, r=50;

  double smx=0, smy=0, smz=0, sz = 0;
  for (int i=0; i<zsol.size(); i++){

    double mx,my,mz;
    if ((i-x0)*(i-x0) + (H-y0)*(H-y0) <= r*r){
      mx = (i-x0)/(double)r;
      my = -(H-y0)/(double)r;
      mz = sqrt(1-mx*mx-my*my);
    }
    else {
      mx = usol[i*NPDE+0];
      my = usol[i*NPDE+1];
      mz = usol[i*NPDE+2];
    }
    int c = vec_to_cal(mx,my,mz);

    ss << ((char)((c>>16)&0xFF))
      << ((char)((c>>8)&0xFF))
      << ((char)(c&0xFF));
  }
  char c=0;
  ss << c << c << c;

  for (int i=0; i<zsol.size(); i++){
    double nx = usol[i*NPDE+3];
    double ny = usol[i*NPDE+4];
    double nz = usol[i*NPDE+5];
    int c = vec_to_cal(nx,ny,nz);
    ss << ((char)((c>>16)&0xFF))
       << ((char)((c>>8)&0xFF))
       << ((char)(c&0xFF));
  }
  H++;
  ss.seekp(width_pos);
  ss << H;
}

void
pnm_writer::hline(){

  // open file for appending; write to its end
  std::ofstream ss;
  std::ios_base::openmode flags =
     (H==0)? std::fstream::out | std::fstream::trunc:
             std::fstream::out | std::fstream::in | std::fstream::ate;
  ss.open (fname, flags);

  for (int i=0; i<W; i++){
    char c = (i/10)%2;
    ss << c << c << c;
  }
  H++;
  ss.seekp(width_pos);
  ss << H;
}

