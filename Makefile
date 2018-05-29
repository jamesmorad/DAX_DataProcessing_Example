#---------------------------------------------------
all:
	g++ -Wall -O2 `root-config --cflags` GenerateRQs.cpp -o GenerateRQs `root-config --libs`

clean:
	rm -rf GenerateRQs


