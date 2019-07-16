CFLAGS    += -g -I../include -I.
LDFLAGS   += -g -L../lib -lcmat -lcunit

CSRC      := main.c \
             test_add.c \
             test_sub.c \
             test_mul.c \
             test_product.c \
             test_det.c \
             test_dot.c \
             test_inverse.c \
             test_transpose.c

OBJS      := $(patsubst %.c,%.o, $(CSRC))

TARGET    := test_all

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

test_add.o : test_add.c test_add.h
test_sub.o : test_sub.c test_sub.h
test_mul.o : test_mul.c test_mul.h
test_product.o : test_product.c test_product.h
test_transpose.o : test_transpose.c test_transpose.h
test_det.o : test_det.c test_det.h
test_dot.o : test_dot.c test_dot.h
test_inverse.o : test_inverse.c test_inverse.h

test: $(TARGET)
	./$(TARGET)

test_add.h:
	ruby make_add_data.rb > $@

test_sub.h:
	ruby make_sub_data.rb > $@

test_mul.h:
	ruby make_mul_data.rb > $@

test_product.h:
	ruby make_product_data.rb > $@

test_transpose.h:
	ruby make_transpose_data.rb > $@

test_det.h:
	ruby make_det_data.rb > $@

test_dot.h:
	ruby make_dot_data.rb > $@

test_inverse.h:
	ruby make_inverse_data.rb > $@

clean:
	rm -f $(TARGET) $(OBJS)

clean-data:
	rm test_*.h
