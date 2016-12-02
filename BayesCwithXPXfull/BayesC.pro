TEMPLATE = app
OBJECTS_DIR   = ./objectDir
CONFIG       = release
DEPENDPATH   += .
INCLUDEPATH  +=  . /Users/erxingfangshui/Dropbox/eigen3/Eigen \
                   /Users/erxingfangshui/Dropbox/matvec 
LIBS +=            /Users/erxingfangshui/Dropbox/matvec/lib/libmatvec.a

# Input
HEADERS += 
SOURCES += mainZPZ.cpp

QMAKE_CXXFLAGS += -O3 -arch x86_64
QMAKE_CXXFLAGS += -funsafe-math-optimizations

QMAKE_CXX =  mpic++
QMAKE_LINK = mpic++	  
