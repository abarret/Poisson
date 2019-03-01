SAMRAI	      =	${HOME}/sfw/samrai/3.12.0/SAMRAI#/not_backed_up/shared/samrai/master/SAMRAI
OBJECT        = ${HOME}/sfw/samrai/3.12.0/linux-dbg#/not_backed_up/shared/samrai/master/linux-dbg

include $(OBJECT)/config/Makefile.config

CPPFLAGS_EXTRA= -DDISPLAY
#LDFLAGS_EXTRA= --link_command_prefix "/usr/local/bin/purify4.5 -best-effort"

OBJS      = main.o CartSideDoubleRefineStrategy.o Poisson.o
main2d:		$(OBJS) $(LIBSAMRAI)
		$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) \
		$(LIBSAMRAI) $(LDLIBS) -o main2d

clean:
		$(RM) main2d core
		$(RM) *.o *.ii *.int.c stamp-[23]d 
		$(RM) -r ti_files ii_files
		$(RM) *log*

include Makefile.depend
