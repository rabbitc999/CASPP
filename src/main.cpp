// SPDX-License-Identifier: MIT
// Copyright (c) 2026 小兔糖糖OvO

#include "cas/engine.h"

#include <cstdio>
#include <cstring>
#include <string>

#ifdef _WIN32
#include <io.h>
#define CAS_ISATTY _isatty
#define CAS_FILENO _fileno
#else
#include <unistd.h>
#define CAS_ISATTY isatty
#define CAS_FILENO fileno
#endif

namespace {

/** @brief 版本信息。 */
const char *const kVersion = "1.0.0";

/**
 * @brief 去掉结果文本里的内部标记。
 *
 * 矩阵结果在 `Result::text` 开头带一个 `\x1F` 标记，供嵌套调用（如 `det(inv(A))`）
 * 识别"上一步的结果是个矩阵"；打印给用户之前要剥掉。
 */
std::string displayText(const std::string &text) {
    std::size_t begin = 0;
    while (begin < text.size() && text[begin] == '\x1F') ++begin;
    return text.substr(begin);
}

/** @brief 打印横幅。 */
void printBanner() {
    std::printf("caspp %s —— 代数符号计算器（C++17）\n", kVersion);
    std::printf("输入 help 查看命令，quit 退出。\n\n");
}

/** @brief 打印单条结果。 */
void printResult(const Cas::Result &result) {
    if (result.ok) {
        std::printf("%s\n", displayText(result.text).c_str());
        std::printf("  [LaTeX] %s\n", result.latex.c_str());
    } else {
        std::printf("错误：%s\n", result.message.c_str());
        if (result.position > 0) std::printf("      位置 %zu\n", result.position);
    }
    std::printf("  [%.2f ms]\n", result.elapsedMs);
}

/** @brief 读取一行（去掉换行）。返回 false 表示到达 EOF。 */
bool readLine(const char *prompt, std::string *line) {
    std::fputs(prompt, stdout);
    std::fflush(stdout);

    std::string buffer;
    int c = 0;
    while ((c = std::fgetc(stdin)) != EOF) {
        if (c == '\n') break;
        if (c != '\r') buffer.push_back(static_cast<char>(c));
    }
    if (c == EOF && buffer.empty()) return false;
    *line = buffer;
    return true;
}

/** @brief 处理一条 REPL 命令；返回 false 表示退出。 */
bool handleCommand(Cas::Engine *engine, const std::string &line) {
    std::string trimmed = line;
    // 去首尾空白
    const auto begin = trimmed.find_first_not_of(" \t");
    if (begin == std::string::npos) return true;
    const auto end = trimmed.find_last_not_of(" \t");
    trimmed = trimmed.substr(begin, end - begin + 1);

    if (trimmed == "quit" || trimmed == "exit" || trimmed == "q") return false;
    if (trimmed == "help" || trimmed == "?") {
        std::printf("%s\n", Cas::Engine::helpText().c_str());
        return true;
    }
    if (trimmed == "unicode on") {
        engine->setUnicodeOutput(true);
        std::printf("输出将使用 Unicode 字符。\n");
        return true;
    }
    if (trimmed == "unicode off") {
        engine->setUnicodeOutput(false);
        std::printf("输出使用 ASCII 字符。\n");
        return true;
    }
    if (trimmed == "mem" || trimmed == "memory") {
        const Cas::MemoryStats stats = engine->memoryStats();
        std::printf("池容量 %zu 字节，活动 %zu 字节（%zu 块），峰值 %zu 字节，累计分配 %zu 次\n",
                    stats.capacity, stats.liveBytes, stats.liveBlocks, stats.peakBytes, stats.totalBlocks);
        return true;
    }
    if (trimmed == "%") {
        std::printf("%s\n", displayText(engine->lastResult()).c_str());
        return true;
    }

    printResult(engine->evaluate(trimmed));
    return true;
}

/** @brief 非交互模式：逐行读 stdin 或执行参数里的表达式。 */
int runBatch(Cas::Engine *engine, const std::vector<std::string> &expressions) {
    if (!expressions.empty()) {
        for (const std::string &expression : expressions) {
            std::printf("> %s\n", expression.c_str());
            printResult(engine->evaluate(expression));
        }
        return 0;
    }

    std::string line;
    while (readLine("", &line)) {
        if (!line.empty() && line[0] == '#') continue;
        std::printf("%s\n", line.c_str());
        if (!handleCommand(engine, line)) break;
    }
    return 0;
}

} // namespace

int main(int argc, char **argv) {
    Cas::Engine engine;

    std::vector<std::string> expressions;
    bool interactive = true;

    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        if (argument == "-h" || argument == "--help") {
            std::printf("用法：caspp [选项] [表达式...]\n\n");
            std::printf("  -c, --command <expr>  求值并退出（可重复）\n");
            std::printf("  -q, --quiet           不打印横幅\n");
            std::printf("  -u, --unicode         使用 Unicode 输出（π √ ·）\n");
            std::printf("  -v, --version         显示版本\n");
            std::printf("  -h, --help            显示本帮助\n\n");
            std::printf("%s\n", Cas::Engine::helpText().c_str());
            return 0;
        }
        if (argument == "-v" || argument == "--version") {
            std::printf("caspp %s\n", kVersion);
            return 0;
        }
        if (argument == "-u" || argument == "--unicode") {
            engine.setUnicodeOutput(true);
            continue;
        }
        if (argument == "-q" || argument == "--quiet") {
            interactive = false;
            continue;
        }
        if (argument == "-c" || argument == "--command") {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "错误：%s 需要一个表达式参数\n", argument.c_str());
                return 2;
            }
            expressions.push_back(argv[++i]);
            interactive = false;
            continue;
        }
        expressions.push_back(argument);
        interactive = false;
    }

    if (!expressions.empty()) return runBatch(&engine, expressions);

    // 无参数：stdin 是终端则进 REPL，否则按批处理读管道
    if (CAS_ISATTY(CAS_FILENO(stdin))) {
        if (interactive) printBanner();
        std::string line;
        while (readLine("cas> ", &line)) {
            if (!handleCommand(&engine, line)) break;
        }
        std::printf("\n");
        return 0;
    }
    return runBatch(&engine, expressions);
}
