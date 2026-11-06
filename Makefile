CC      := gcc
CFLAGS  := -g -O2 -I src/include -I src/SDL3_gfx -MMD -MP 
LDFLAGS := -L src/lib
LDLIBS  := -lSDL3 -lSDL3_ttf -lm

SRCS := GravSim.c $(wildcard ./src/*.c)  $(wildcard ./src/SDL3_gfx/*.c)

OBJS := $(addprefix build/,$(notdir $(SRCS:.c=.o)))
DEPS := $(OBJS:.o=.d)
vpath %.c . ./src ./src/SDL3_gfx

ifeq ($(OS),Windows_NT)
    TARGET := GravSim.exe
    MKDIR = if not exist "$(subst /,\,$(1))" mkdir "$(subst /,\,$(1))"
    CLEAN = if exist build rmdir /s /q build & if exist $(TARGET) del /q $(TARGET)
else
    TARGET := GravSim
    MKDIR = mkdir -p $(1)
    CLEAN = rm -rf build $(TARGET)
endif

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

build/%.o: %.c
	@$(call MKDIR,build)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(CLEAN)

.PHONY: all clean

-include $(DEPS)