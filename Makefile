# ---------------------------------------------------------------------------
# caspp —— C++17 代数符号计算器
#
#   make              构建 build/caspp（可执行）与全部测试
#   make test         运行全部测试（等价于 make check）
#   make run          启动 REPL
#   make bench        跑大整数与表达式性能基准
#   make clean        清理 build/
#
# 目录约定：头文件 include/cas/**，实现 src/**（对应目录结构），测试 tests/**。
# 产物全部落在 build/，源码树保持干净。
# ---------------------------------------------------------------------------

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic
CPPFLAGS += -Iinclude
LDFLAGS  ?=
LDLIBS   ?=

BUILD_DIR := build

# 库源文件（不含 main.cpp）
LIB_SOURCES := \
	$(filter-out src/main.cpp,$(wildcard src/*.cpp)) \
	$(wildcard src/detail/*.cpp)
LIB_OBJECTS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(LIB_SOURCES))

# 可执行与测试
MAIN_OBJECT := $(BUILD_DIR)/src/main.o
TARGET      := $(BUILD_DIR)/caspp

# 测试只链接库对象 + 各自的 main，不能带上 src/main.o
# 注意用递归赋值：LIB_OBJECTS 在其后才定义
TEST_OBJECTS = $(LIB_OBJECTS)

TEST_SOURCES := $(wildcard tests/*.cpp)
TEST_TARGETS := $(patsubst tests/%.cpp,$(BUILD_DIR)/tests/%,$(TEST_SOURCES))

.PHONY: all test check run bench clean rebuild format

all: $(TARGET) $(TEST_TARGETS)

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(TARGET): $(LIB_OBJECTS) $(MAIN_OBJECT)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)
	@echo "  -> $@"

$(BUILD_DIR)/tests/%: tests/%.cpp $(TEST_OBJECTS)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(TEST_OBJECTS) -o $@ $(LDFLAGS) $(LDLIBS)

# 头文件依赖：简单带过（改头文件时整体重编，避免漏依赖）
$(LIB_OBJECTS) $(MAIN_OBJECT): $(wildcard include/cas/*.h include/cas/detail/*.h)

test check: all
	@echo "=== run all tests ==="
	@status=0; \
	for t in $(TEST_TARGETS); do \
		echo "--- $$t"; \
		$$t || status=1; \
	done; \
	if [ $$status -eq 0 ]; then echo "=== all passed ==="; else echo "=== FAILURES ==="; fi; \
	exit $$status

run: $(TARGET)
	@./$(TARGET)

bench: $(TARGET)
	@./$(BUILD_DIR)/tests/bigIntFuzz >/dev/null 2>&1 || true
	@echo "=== 基准 ==="
	@./$(TARGET) -q -c "factor(x^12-1)" -c "expand((x+y+z)^8)" -c "integrate(x^3*sin(x), x)" \
		-c "groebner([x^2+y^2-1, x-y], [x, y])" -c "solve(x^4-5x^2+4=0, x)"

clean:
	rm -rf $(BUILD_DIR)

rebuild: clean all
