CFLAGS    += -O3 -g -I../include -I.
LDFLAGS   += -g -L../lib -lcmat -lcunit

CSRC      := main.c \
             test_new.c \
             test_clone.c \
             test_destroy.c \
             test_append.c \
             test_add.c \
             test_sub.c \
             test_mul.c \
             test_product.c \
             test_det.c \
             test_dot.c \
             test_inverse.c \
             test_transpose.c \
             test_lu_decomp.c \
	     test_permute_row.c \
	     test_permute_column.c

OBJS      := $(patsubst %.c,%.o, $(CSRC))

HEADERS   := test_add.h \
             test_sub.h \
             test_mul.h \
             test_product.h \
             test_transpose.h \
             test_det.h \
             test_dot.h \
             test_inverse.h \
             test_lu_decomp.h

TARGET    := test_all

all: $(TARGET)

$(TARGET): $(OBJS) ../lib/libcmat.a
	$(CC) -o $@ $^ $(LDFLAGS)

benchmark: benchmark.c ../lib/libcmat.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

tt: tt.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)


test_new.o : test_new.c
test_clone.o : test_clone.c
test_destroy.o : test_destroy.c
test_append.o : test_append.c
test_add.o : test_add.c test_add.h
test_sub.o : test_sub.c test_sub.h
test_mul.o : test_mul.c test_mul.h
test_product.o : test_product.c test_product.h
test_transpose.o : test_transpose.c test_transpose.h
test_det.o : test_det.c test_det.h
test_dot.o : test_dot.c test_dot.h
test_inverse.o : test_inverse.c test_inverse.h
test_lu_decomp.o : test_lu_decomp.c test_lu_decomp.h
test_permute_row.o: test_permute_row.c
test_permute_column.o: test_permute_column.c

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

test_lu_decomp.h:
	ruby make_lu_decomp_data.rb > $@

data: $(HEADERS)

clean:
	rm -f $(TARGET) $(OBJS) benchmark

clean-data:
	rm test_*.h
