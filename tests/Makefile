CFLAGS    += -I../include -I.
LDFLAGS   += -L../lib -lcmat -lcunit

CSRC      := main.c \
             test_det.c \
             test_inverse.c \
             test_mul.c \
             test_transpose.c

OBJS      := $(patsubst %.c,%.o, $(CSRC))

TARGET    := test_all

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

main.c : test.h

test: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) $(OBJS)
