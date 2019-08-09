CFLAGS    += -O0 -g -I./include -DDEBUG
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

clean:
	make -C tests $@
	rm -rf $(dir $(TARGET)) $(OBJS)
