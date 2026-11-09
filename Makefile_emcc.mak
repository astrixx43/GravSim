CC      := emcc
BUILD   := build-web
PORTS   := -sUSE_SDL=3 -sUSE_SDL_TTF=3
CFLAGS  := -O2 -msimd128 -msse -I ./src/include -I ./src/SDL3_gfx $(PORTS) -MMD -MP
LDFLAGS := $(PORTS) -sALLOW_MEMORY_GROWTH=1 --preload-file ./src/fonts/dejavu-sans-condensed-webfont.ttf
LDLIBS  := -lm

SRCS := emcc_project/main_emcc.c $(wildcard ./src/*.c) $(wildcard ./src/SDL3_gfx/*.c)
OBJS := $(addprefix $(BUILD)/,$(notdir $(SRCS:.c=.o)))
DEPS := $(OBJS:.o=.d)
vpath %.c ./emcc_project ./src ./src/SDL3_gfx

TARGET := ./emcc_project/GravSim.html

ifeq ($(OS),Windows_NT)
    MKDIR = if not exist "$(subst /,\,$(1))" mkdir "$(subst /,\,$(1))"
    CLEAN = if exist $(BUILD) rmdir /s /q $(BUILD) & del /q GravSim.html GravSim.js GravSim.wasm GravSim.data 2>nul
else
    MKDIR = mkdir -p $(1)
    CLEAN = rm -rf $(BUILD) GravSim.html GravSim.js GravSim.wasm GravSim.data
endif

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(BUILD)/%.o: %.c
	@$(call MKDIR,$(BUILD))
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(CLEAN)

.PHONY: all clean
-include $(DEPS)