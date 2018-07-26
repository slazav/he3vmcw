FFLAGS= -Werror -Wconversion\
  -Wintrinsic-shadow -Wline-truncation\
  -Waliasing  -Wampersand -Warray-bounds -Wcharacter-truncation\
  -Wline-truncation -Wintrinsics-std -Wsurprising -Wno-tabs -Wunderflow\
  -Wno-unused-parameter -Wno-align-commons

TARGETS=vmcw vft

all: $(TARGETS)
clean:
	rm -f $(TARGETS)

vmcw: vmcw.f libs/pde_dp.f he3_funcs.f\
      vmcw_mesh.f vmcw_mon.f vmcw_state.f vmcw_pdecol.f\
      vmcw_cmd.f vmcw_func.f
vft: vft.f

