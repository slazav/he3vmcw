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

/****************************************************************/

pnm_writer_t::pnm_writer_t(
    const std::string & fname_,
    const pdecol_solver *s,
    const writer_type_t type_):
  H(0),W(0),fname(fname_),solver(s),
  br(0), bx0(0), by0(0), hline_n(0), type(type_){
}

void
pnm_writer_t::write(const std::vector<double> & zsol,
                    const std::vector<double> & usol){
  if (!solver) return;

  int NPDE = solver->get_npde();

  // open file for appending; write to its end
  std::ofstream ss;
  std::ios_base::openmode flags =
     (H==0)? std::fstream::out | std::fstream::trunc :
             std::fstream::out | std::fstream::in | std::fstream::ate;
  ss.open (fname, flags);

  int NPTS = zsol.size();
  int NDER = usol.size()/NPTS/NPDE;
  if (NDER<1) return;

  if (type==WRITER_PNM) {
    // if it is the first line
    if (H==0){
      ss.seekp(0);
      W = NPTS*2+1;
      ss << "P6\n" << W << " ";
      width_pos = ss.tellp();
      ss << "                  \n255\n";
    }

    char ck=0;
    for (; hline_n>0; hline_n--){
      for (int i=0; i<W; i++) ss << ck << ck << ck;
      H++;
    }

    double smx=0, smy=0, smz=0, sz = 0;
    for (int i=0; i<NPTS; i++){

      double mx,my,mz;
      if (br && (i-bx0)*(i-bx0) + (H-by0)*(H-by0) <= br*br){
        mx = (i-bx0)/(double)br;
        my = -(H-by0)/(double)br;
        mz = sqrt(1-mx*mx-my*my);
      }
      else {
        mx = usol[i*NPDE+0];
        my = usol[i*NPDE+1];
        mz = usol[i*NPDE+2];
      }
      int col = vec_to_cal(mx,my,mz);

      ss << ((char)((col>>16)&0xFF))
         << ((char)((col>>8)&0xFF))
         << ((char)(col&0xFF));
    }
    ss << ck << ck << ck;

    for (int i=0; i<NPTS; i++){
      double nx = usol[i*NPDE+3];
      double ny = usol[i*NPDE+4];
      double nz = usol[i*NPDE+5];
      int col = vec_to_cal(nx,ny,nz);
      ss << ((char)((col>>16)&0xFF))
         << ((char)((col>>8)&0xFF))
         << ((char)(col&0xFF));
    }
    H++;
    ss.seekp(width_pos);
    ss << H;
    return;
  }

  if (type==WRITER_TXT){

    // print header
    if (H==0) ss << "# n  M\n";

    // print blank lines
    for (; hline_n>0; hline_n--){ ss << "\n"; }

    double smx=0, smy=0, smz=0, sz = 0;
    for (int i=0; i<NPTS; i++){

      double z  = zsol[i];
      double mx = usol[i*NPDE+0];
      double my = usol[i*NPDE+1];
      double mz = usol[i*NPDE+2];

      double nx = usol[i*NPDE+3];
      double ny = usol[i*NPDE+4];
      double nz = usol[i*NPDE+5];
    }
    ss << "\n";
    return;
  }

}

void
pnm_writer_t::hline(){ hline_n++; }

void
pnm_writer_t::legend(int r, int x0){
  br = r; bx0=x0; by0 = H+r;
}

/****************************************************************/

bool
pnm_writers_t::add(const std::string & fname,
                   const pdecol_solver *solver,
                   const writer_type_t type){
   iterator w = find(fname);
   if (w != end()) return false;
   insert(std::make_pair(fname, pnm_writer_t(fname, solver, type)));
  return true;
}

bool
pnm_writers_t::del(const std::string & fname){
   iterator w = find(fname);
   if (w == end()) return false;
   erase(w);
   return true;
}

bool
pnm_writers_t::legend(const std::string & fname, int r, int x0){
  iterator w = find(fname);
  if (w == end()) return false;
  w->second.legend(r,x0);
  return true;
}

bool
pnm_writers_t::hline(const std::string & fname){
  iterator w = find(fname);
  if (w == end()) return false;
  w->second.hline();
  return true;
}


void
pnm_writers_t::write(const std::vector<double> zsol,
           const std::vector<double> usol){
  for (iterator i=begin(); i!=end(); i++)
    i->second.write(zsol, usol);
}



