#ifndef PNM_WRITER_H
#define PNM_WRITER_H

#include <vector>
#include <string>

class pnm_writer_t{
  public:
    pnm_writer_t(const std::string & fname);

    // write line of data
    void write(const std::vector<double> zsol,
               const std::vector<double> usol, int NPDE);

    // draw a horisontal line (to mark some position)
    void hline();

    // start drawing a color legend with radius r and x position x0
    void legend(int r=50, int x0=50);

  private:
    std::string fname;
    int width_pos;
    int W,H;
    int br, bx0, by0;
    int hline_n;
};

#include <map>

// Container for pnm_writers. Key is file name.
class pnm_writers_t : public std::map<std::string, pnm_writer_t>{
  public:

    // Add a pnm_writer for file fname, return false if
    // it already exists.
    bool add(const std::string & fname);

    // Delete a pnm_writer for file fname. Return false if
    // pnm_writer did not exist, true if it has been deleted
    bool del(const std::string & fname);

    // Run legend command for fname pnm_writer. Return false if
    // pnm_writer do not exist, true if command has been run.
    bool legend(const std::string & fname, int r=50, int x0=50);

    // Run hline command for fname pnm_writer. Return false if
    // pnm_writer do not exist, true if command has been run.
    bool hline(const std::string & fname);

    // write line of data by all writers
    void write(const std::vector<double> zsol,
               const std::vector<double> usol, int NPDE);

};


#endif