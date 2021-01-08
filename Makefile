CXXC = g++

STD = -std=c++17

LDFLAGS = -lm

OPT = -D NDEBUG -ffast-math -O3 -march=native

DBG = -D DEBUG -Og -g -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy \
                      -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op          \
                      -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept             \
                      -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow     \
                      -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel                \
                      -Wstrict-overflow=5 -Wswitch-default -Wundef -Werror -Wno-unused     \

SRC =                        \
./source/main.cpp            \
./source/paut.cpp            \
./source/readData.cpp        \
./source/writeLawFile.cpp    \
./source/solveCylSnell.cpp   \
./source/refractedAngle.cpp  \
./source/sourcePoints.cpp    \
./source/exitPoints.cpp      \
./source/focalAbberation.cpp \
./source/root_rc.cpp         \
./source/timeDelays.cpp      \
./source/reachingTimes.cpp   \
./source/optimizeFocus.cpp   \
./source/local_min_rc.cpp    \
./source/gDelay.cpp          \

INCL = -I ./source/headers

run: $(SRC)
	$(CXXC) $(STD) $(OPT) $(INCL) $(SRC) $(LDFLAGS) -o run

dbg: $(SRC)
	$(CXXC) $(STD) $(DBG) $(INCL) $(SRC) $(LDFLAGS) -o dbg

clear:
	rm -f run dbg
