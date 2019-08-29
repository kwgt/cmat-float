CFLAGS    += -g -I./include -DDEBUG

CFLAGS    += -O3 -ftree-vectorize
CFLAGS    += -DENABLE_NEON
CFLAGS    += -fopenmp

LDFLAGS   += -g -L./lib -lcunit

CSRC      := src/cmat.c
OBJS      := $(patsubst %.c,%.o, $(CSRC))

TARGET    := lib/libcmat.a

all: $(TARGET)

$(TARGET): $(OBJS)
	@if [ ! -e $(dir $(TARGET)) ];then mkdir -p $(dir $(TARGET));fi
	ar rcs $@ $^
	ranlib $@

src/cmat.c: include/cmat.h


test:
	make -C tests

benchmark:
	make -C tests benchmark

clean:
	rm -rf $(dir $(TARGET)) $(OBJS)

clean-test:
	make -C tests $@
