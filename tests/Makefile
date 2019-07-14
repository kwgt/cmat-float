CFLAGS    += -g -I../include -I.
LDFLAGS   += -g -L../lib -lcmat -lcunit

CSRC      := main.c \
             test_add.c \
             test_sub.c \
             test_mul.c \
             test_det.c \
             test_inverse.c \
             test_transpose.c

OBJS      := $(patsubst %.c,%.o, $(CSRC))

TARGET    := test_all

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

test: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) $(OBJS)
