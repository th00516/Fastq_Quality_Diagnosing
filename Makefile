CXX      = g++
CXXFLAGS = -lboost_iostreams
EXEC     = fastqQualDiagnosing_lite

$(EXEC): main.cpp fastqQualDiagnosing_lite.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

clean:
	rm $(EXEC)
