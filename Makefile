all: ThiessenInterp

ThiessenInterp: ConsoleApplication2.o OGRFile.o st2ws.o
	 g++ -o ./Debug/ThiessenInterp -L/usr/local/lib  ConsoleApplication2.o OGRFile.o st2ws.o -lgeos -lgdal  `wx-config --libs all` -lpthread

ConsoleApplication2.o:ConsoleApplication2.cpp
	 g++ -c -std=c++11 -I/usr/local/include/wx-3.1 -I/usr/local/include -g `wx-config --cxxflags --libs std` ConsoleApplication2.cpp

OGRFile.o:OGRFile.cpp
	 g++ -c -std=c++11 -g OGRFile.cpp
st2ws.o:st2ws.cpp
	 g++ -c -std=c++11 -g st2ws.cpp
     
clean:
	 rm ConsoleApplication2.o OGRFile.o st2ws.o ./Debug/ThiessenInterp
