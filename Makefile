FFLAGS= -Werror -Wconversion\
  -Wintrinsic-shadow -Wline-truncation\
  -Waliasing  -Wampersand -Warray-bounds -Wcharacter-truncation\
  -Wline-truncation -Wintrinsics-std -Wsurprising -Wno-tabs -Wunderflow\
  -Wno-unused-parameter -Wno-align-commons

all: vmcw vmcw_st he3
vmcw: vmcw.f libs/pde_dp.f he3_funcs.f
he3: he3.f he3_funcs.f
