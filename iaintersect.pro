######################################################################
# pro file for complexity plot Wed Oct 19 14:59:16 2011
######################################################################

TEMPLATE = app
TARGET   = iaintersect

CONFIG   += console warn_on release
CONFIG   -= app_bundle

QT       -= gui

!win32{
DEFINES                  += _APPNAME=\\\"$$TARGET\\\"
OBJECTS_DIR              = GeneratedFiles
UI_DIR                   = GeneratedFiles
MOC_DIR                  = GeneratedFiles
RCC_DIR                  = GeneratedFiles
INCLUDEPATH += /usr/local/include/
}

macx{
QMAKE_CFLAGS_X86_64 += -mmacosx-version-min=10.7
QMAKE_CXXFLAGS_X86_64 = $$QMAKE_CFLAGS_X86_64
}

win32{
DEFINES        += _APPNAME=\"$$TARGET\"
}

DEPENDPATH  +=  . \
                ./src

INCLUDEPATH +=  . \
                ./src

HEADERS     +=  ./src/iaintersect.hpp \
                ./src/Arguments.hpp \
                ./src/main.hpp \
                ./src/config.hpp \
                ./src/peak_reader.hpp \
                ./src/string_tools.hpp

SOURCES     +=  ./src/main.cpp \
                ./src/iaintersect.cpp \
                ./src/Arguments.cpp \
                ./src/peak_reader.cpp

QMAKE_CLEAN += $$TARGET logfile.log *~ *.txt
