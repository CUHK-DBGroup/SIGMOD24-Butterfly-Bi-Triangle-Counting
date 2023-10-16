OBJECTS = app.o bi_graph.o util.o tracker.o exact_alg.o gt_manager.o

EXE = app

FLAGS = -Wall -Wextra -O3 -std=c++17

$(EXE): $(OBJECTS)
	g++ $(OBJECTS) -o $(EXE)

-include $(OBJECTS:.o=.d)

%.o: %.cpp
	g++ -c $(FLAGS) $*.cpp -o $*.o
	g++ -MM $(FLAGS) $*.cpp > $*.d

.PHONY: cleanobj
cleanobj:
	rm -f $(OBJECTS)

.PHONY: clean
clean:
	rm -f $(OBJECTS) *.d $(EXE)

