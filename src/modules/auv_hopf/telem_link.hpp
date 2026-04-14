#pragma once

#include <stdint.h>

// --- 1. 定义一个结构体，用于存放你想提取的数据 ---
struct TelemData {
    float rpm{0.0f};         // 举个例子，假设你传的是 RPM
    // float external_x{0.0f}; // 你可以随便加你要提取的变量
    bool updated{false};     // 数据更新标志位
};

class TelemLink
{
public:
    explicit TelemLink(const char *device_path = "/dev/ttyS1");
    ~TelemLink();

    bool init();
    void update();
    void send_string(const char *str);

    // --- 2. 新增一个提取数据的接口 ---
    // 如果有新数据，返回 true 并把数据拷贝给你，然后清空更新标志
    bool get_latest_data(TelemData &out_data);

private:
    const char *_device_path;
    int _fd{-1};

    char _rx_line[128]{};
    int _rx_len{0};

    // --- 3. 新增一个内部变量存数据 ---
    TelemData _latest_data{};

    void handle_uart_line(const char *line);
};
