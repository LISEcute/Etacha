TARGET = ETACHA4
TEMPLATE = app
QT       += core gui printsupport charts
CONFIG += c++17

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

UI_DIR = ui
MOC_DIR = moc
RCC_DIR = rcc
OBJECTS_DIR = obj


win32-g++ {
DESTDIR = c:/Etacha/_install
}
win32-msvc {
DESTDIR = c:/Etacha/_install_MSVC
}

win32:VERSION = 4.4.18.1 # major.minor.patch.build
else:VERSION = 4.4.18    # major.minor.patch

win32 {
       QMAKE_TARGET_COPYRIGHT = "LISE group at FRIB/MSU"
	QMAKE_TARGET_COMPANY   = "LISE group at FRIB/MSU"
	}

#==================================================
# The following define makes your compiler warn you if you use any
# feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

# Input



SOURCES += \
       L_Loss/L_Ziegler.cpp \
       g_Etacha4/source/e_AtomicShell.cpp \
       g_Etacha4/source/e_Auger4.cpp \
       g_Etacha4/source/e_Donaut4.cpp \
       g_Etacha4/source/e_Etacha4.cpp \
       g_Etacha4/source/e_F4.cpp \
       g_Etacha4/source/e_INTG.cpp \
       g_Etacha4/source/e_INTG2.cpp \
       g_Etacha4/source/e_Pion.cpp \
       g_Etacha4/source/e_SEIK.cpp \
       g_Etacha4/source/e_SecMean4.cpp \
       g_Etacha4/source/e_Sex2.cpp \
       g_Etacha4/source/e_Sexi.cpp \
       g_Etacha4/source/e_Snl.cpp \
       g_Etacha4/source/e_Tceis.cpp \
       g_Etacha4/source/e_senlm4.cpp \
       g_Etacha4/util/gl_Element.cpp \
       g_Etacha4/util/gl_Range.cpp \
       g_Etacha4/util/gl_Util.cpp \
       g_Etacha4/win/e_CSresults.cpp \
       g_Etacha4/win/e_actionsPrintAbout.cpp \
       g_Etacha4/win/e_fileReadWrite.cpp \
       g_Etacha4/win/e_graph.cpp \
       g_Etacha4/win/e_graphEvolution.cpp \
       g_Etacha4/win/e_main.cpp\
       g_Etacha4/win/e_about.cpp \
       g_Etacha4/win/e_mainwindow.cpp \
       g_Etacha4/win/e_winUtil.cpp \
    o_ODE/ode.cpp \
    o_ODE/ode_rkf_util.cpp \
    o_ODE/rkf45.cpp \
    w_Stuff/w_Label_clickable.cpp \
    w_Stuff/win_utilString.cpp

HEADERS  += \
       g_Etacha4/source/e_AtomicShell.h \
       g_Etacha4/source/e_Declare.h \
       g_Etacha4/source/e_Etacha4.h \
       g_Etacha4/source/e_GaussDataBlock.h \
       g_Etacha4/source/e_TargetData.h \
       g_Etacha4/util/string_utils.h \
       g_Etacha4/win/e_CSresults.h \
       g_Etacha4/win/e_Constant.h \
       g_Etacha4/win/e_ftype.h \
	g_Etacha4/win/e_about.h \
       g_Etacha4/win/e_graph.h \
       g_Etacha4/win/e_mainwindow.h \
       g_Etacha4/win/e_myextern.h \
	o_ODE/ode.hpp \
	o_ODE/rkf45.hpp \
	w_Stuff/liseStrcpyOS.h \
	w_Stuff/w_Label_clickable.h

FORMS    += \
       g_Etacha4/win/e_CSresults.ui \
	g_Etacha4/win/e_about.ui \
       g_Etacha4/win/e_graph.ui \
       g_Etacha4/win/e_mainwindow.ui

RESOURCES += \
       g_Etacha4/Icons/e_Icons.qrc \
       lise.qrc

RC_ICONS = g_Etacha4/Icons/etacha.ico

# probably "macx" instead "mac"
#mac {
#DEFINES += __APPLE_
#ICON = ./Icons_macos/etacha.icns
#QMAKE_INFO_PLIST = ./Info.plist
#installFolder.files = _install
#installFolder.path = Contents
#QMAKE_BUNDLE_DATA += installFolder
#}
