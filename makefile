TARGET=a.out

CC=g++
CPPFLAGS=-O3 -std=c++11

$(TARGET): thermostats.cpp
	$(CC) $(CPPFLAGS) $< -o $@

clean:
	rm -f $(TARGET) *.dat *.png
