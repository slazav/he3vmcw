#ifndef PNM_WRITER_H
#define PNM_WRITER_H

#include <iostream>
#include <vector>

class pnm_writer{
  public:
    pnm_writer(const char *fname_);

    // write line of data
    void write(const std::vector<double> zsol,
               const std::vector<double> usol, int NPDE);

    // draw a horisontal line (to mark some position)
    void hline();

    // start drawing a color legend with radius r and x position x0
    void legend(int r=50, int x0=50);

  private:
    const char *fname;
    int width_pos;
    int W,H;
    int br, bx0, by0;
};

#endif