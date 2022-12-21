NAME        = ../../test2/myprogram
SRCS        = ./src/*cpp ./src/potential/*cpp ./src/boundary/*cpp ./src/analysis/*cpp ./src/initialization/*cpp ./src/output/*cpp
OBJS        = $(SRCS:.cpp=.o)
CXX         = g++
CXXFLAGS    = -std=c++11 -fopenmp

.PHONY: all
all: $(NAME)

$(NAME): $(OBJS)
	$(CXX) -O3 $(CXXFLAGS) -o $(NAME) $(OBJS)

.PHONY: clean
clean:
	$(RM) $(OBJS)

.PHONY: fclean
fclean: clean
	$(RM) $(NAME)

.PHONY: re
re: fclean all
