FC      = gfortran
FFLAGS  = -O2 -Wall -ffree-line-length-none
TARGET  = test_vtu
SOURCES = mod_vtu.f90 main.f90
OBJECTS = $(SOURCES:.f90=.o)
OUTPUTS = out_ascii.vtu out_binary.vtu mixed_ascii.vtu mixed_binary.vtu

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# main 依赖模块: 必须先编译 mod_vtu.f90 生成 vtu_writer.mod
main.o: mod_vtu.o

run: $(TARGET)
	./$(TARGET)

# 自检脚本放在上一级的 tools/ 里
check: run
	python3 ../tools/validate_vtu.py .

check-vtk: run
	python3 ../tools/validate_vtu_vtk.py .

clean:
	rm -f $(OBJECTS) $(TARGET) *.mod

distclean: clean
	rm -f $(OUTPUTS)

.PHONY: all run check check-vtk clean distclean
