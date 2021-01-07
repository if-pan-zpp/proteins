# This produces a coverage report in an .html file.
# Just put the cg.f and inputfile in this folder
# (together with anything cg.f might need during runtime) and run make.
# 
# The result will be visible in result/index.html.
# 
# You need to install lcov and check that gcov and gfortran have the same version.
#
# Some lines may have no counts at all, it's because they were optimized out
# (sadly, -O0 breaks the program).

all: 
	gfortran -fprofile-arcs -ftest-coverage -O1 -pg -o cg cg.f
	./cg inputfile
	gcov cg.f
	lcov --gcov-tool gcov --capture --directory . --output-file coverage.info
	genhtml --output-directory result coverage.info

clean:
	rm -rf *.mod *.o *.gcov cg gmon.out

cleanall: clean
	rm -rf result/ *.gcda *.gcno coverage.info
