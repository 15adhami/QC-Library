all:
	${CXX} -std=c++11 main.cpp qc.cpp
clean:
	rm -f a.out