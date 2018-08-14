#ifndef PNM_WRITER_H
#define PNM_WRITER_H

#include <iostream>
#include <vector>

class pnm_writer{
  public:
    pnm_writer(std::ostream & ss_);
    void write(const std::vector<double> zsol,
               const std::vector<double> usol, int NPDE);
  private:
    std::ostream & ss;
    int width_pos;
    int line;
};

#endif