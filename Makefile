#
# Makefile for autoconf tutorial
#

INCLUDE = -I$(HOME)/Software/include
WXCONFIG = `wx-config --cxxflags`
FRAMEWORK = -framework OpenGL
FRAMEWORK += -framework GLUT
FRAMEWORK += -framework Accelerate

CXX = g++
DEFS = -O3 -DNDEBUG $(INCLUDE) $(WXCONFIG)
LBFGSLIB = -llbfgs
WXLIBS = `wx-config --libs --gl-libs`
LIBS = $(FRAMEWORK) -lpthread -llapack -lblas $(LBFGSLIB) $(WXLIBS)

SRCS = DrawPanel.cpp SubPanel.cpp MainFrame.cpp AppearanceWindow.cpp Main.cpp calcLayout.cpp constraintsolver3d.cpp constraintsolver2d.cpp
OBJS = $(SRCS:.cpp=.o)
PROG = agi3d

all: $(PROG)

$(PROG): $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LIBS)

.cpp.o:
	$(CXX) $(DEFS) -c -o $@ $<

clean:
	rm -f $(OBJS)

veryclean: clean
	rm $(PROG)
