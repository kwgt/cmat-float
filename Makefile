CFLAGS    += -O0 -g -I./include
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

tests/test_all: $(TARGET)
	make -C tests $(@F)

test: tests/test_all
	./$<

clean:
	make -C tests $@
	rm -rf $(dir $(TARGET)) $(OBJS)
