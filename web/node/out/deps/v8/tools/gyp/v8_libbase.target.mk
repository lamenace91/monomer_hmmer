# This file is generated by gyp; do not edit.

TOOLSET := target
TARGET := v8_libbase
DEFS_Debug := \
	'-DV8_TARGET_ARCH_IA32' \
	'-DENABLE_DISASSEMBLER' \
	'-DV8_ENABLE_CHECKS' \
	'-DOBJECT_PRINT' \
	'-DVERIFY_HEAP' \
	'-DDEBUG' \
	'-DENABLE_EXTRA_CHECKS' \
	'-DENABLE_HANDLE_ZAPPING' \
	'-D_DEBUG'

# Flags passed to all source files.
CFLAGS_Debug := \
	-pthread \
	-Wall \
	-Wextra \
	-Wno-unused-parameter \
	-m32 \
	-fno-strict-aliasing \
	-m32 \
	-Woverloaded-virtual \
	 \
	-fdata-sections \
	-ffunction-sections \
	-g

# Flags passed to only C files.
CFLAGS_C_Debug :=

# Flags passed to only C++ files.
CFLAGS_CC_Debug := \
	-fno-rtti \
	-fno-exceptions

INCS_Debug := \
	-I$(srcdir)/deps/v8

DEFS_Release := \
	'-DV8_TARGET_ARCH_IA32' \
	'-DENABLE_DISASSEMBLER'

# Flags passed to all source files.
CFLAGS_Release := \
	-pthread \
	-Wall \
	-Wextra \
	-Wno-unused-parameter \
	-m32 \
	-fno-strict-aliasing \
	-m32 \
	-O3 \
	-ffunction-sections \
	-fdata-sections \
	-fno-tree-vrp \
	-fno-omit-frame-pointer \
	-fdata-sections \
	-ffunction-sections \
	 \
	-O3

# Flags passed to only C files.
CFLAGS_C_Release :=

# Flags passed to only C++ files.
CFLAGS_CC_Release := \
	-fno-rtti \
	-fno-exceptions

INCS_Release := \
	-I$(srcdir)/deps/v8

OBJS := \
	$(obj).target/$(TARGET)/deps/v8/src/base/atomicops_internals_x86_gcc.o \
	$(obj).target/$(TARGET)/deps/v8/src/base/cpu.o \
	$(obj).target/$(TARGET)/deps/v8/src/base/logging.o \
	$(obj).target/$(TARGET)/deps/v8/src/base/once.o \
	$(obj).target/$(TARGET)/deps/v8/src/base/platform/time.o \
	$(obj).target/$(TARGET)/deps/v8/src/base/platform/condition-variable.o \
	$(obj).target/$(TARGET)/deps/v8/src/base/platform/mutex.o \
	$(obj).target/$(TARGET)/deps/v8/src/base/platform/semaphore.o \
	$(obj).target/$(TARGET)/deps/v8/src/base/utils/random-number-generator.o \
	$(obj).target/$(TARGET)/deps/v8/src/base/platform/platform-linux.o \
	$(obj).target/$(TARGET)/deps/v8/src/base/platform/platform-posix.o

# Add to the list of files we specially track dependencies for.
all_deps += $(OBJS)

# CFLAGS et al overrides must be target-local.
# See "Target-specific Variable Values" in the GNU Make manual.
$(OBJS): TOOLSET := $(TOOLSET)
$(OBJS): GYP_CFLAGS := $(DEFS_$(BUILDTYPE)) $(INCS_$(BUILDTYPE))  $(CFLAGS_$(BUILDTYPE)) $(CFLAGS_C_$(BUILDTYPE))
$(OBJS): GYP_CXXFLAGS := $(DEFS_$(BUILDTYPE)) $(INCS_$(BUILDTYPE))  $(CFLAGS_$(BUILDTYPE)) $(CFLAGS_CC_$(BUILDTYPE))

# Suffix rules, putting all outputs into $(obj).

$(obj).$(TOOLSET)/$(TARGET)/%.o: $(srcdir)/%.cc FORCE_DO_CMD
	@$(call do_cmd,cxx,1)

# Try building from generated source, too.

$(obj).$(TOOLSET)/$(TARGET)/%.o: $(obj).$(TOOLSET)/%.cc FORCE_DO_CMD
	@$(call do_cmd,cxx,1)

$(obj).$(TOOLSET)/$(TARGET)/%.o: $(obj)/%.cc FORCE_DO_CMD
	@$(call do_cmd,cxx,1)

# End of this set of suffix rules
### Rules for final target.
LDFLAGS_Debug := \
	-pthread \
	-rdynamic \
	-m32 \
	-m32

LDFLAGS_Release := \
	-pthread \
	-rdynamic \
	-m32 \
	-m32

LIBS :=

$(obj).target/deps/v8/tools/gyp/libv8_libbase.a: GYP_LDFLAGS := $(LDFLAGS_$(BUILDTYPE))
$(obj).target/deps/v8/tools/gyp/libv8_libbase.a: LIBS := $(LIBS)
$(obj).target/deps/v8/tools/gyp/libv8_libbase.a: TOOLSET := $(TOOLSET)
$(obj).target/deps/v8/tools/gyp/libv8_libbase.a: $(OBJS) FORCE_DO_CMD
	$(call do_cmd,alink)

all_deps += $(obj).target/deps/v8/tools/gyp/libv8_libbase.a
# Add target alias
.PHONY: v8_libbase
v8_libbase: $(obj).target/deps/v8/tools/gyp/libv8_libbase.a

# Add target alias to "all" target.
.PHONY: all
all: v8_libbase

# Add target alias
.PHONY: v8_libbase
v8_libbase: $(builddir)/libv8_libbase.a

# Copy this to the static library output path.
$(builddir)/libv8_libbase.a: TOOLSET := $(TOOLSET)
$(builddir)/libv8_libbase.a: $(obj).target/deps/v8/tools/gyp/libv8_libbase.a FORCE_DO_CMD
	$(call do_cmd,copy)

all_deps += $(builddir)/libv8_libbase.a
# Short alias for building this static library.
.PHONY: libv8_libbase.a
libv8_libbase.a: $(obj).target/deps/v8/tools/gyp/libv8_libbase.a $(builddir)/libv8_libbase.a

# Add static library to "all" target.
.PHONY: all
all: $(builddir)/libv8_libbase.a
