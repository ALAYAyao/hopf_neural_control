#include "telem_link.hpp"

#include <px4_platform_common/log.h>
#include <fcntl.h>
#include <termios.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>

TelemLink::TelemLink(const char *device_path) :
    _device_path(device_path)
{
}

TelemLink::~TelemLink()
{
    if (_fd >= 0) {
        ::close(_fd);
        _fd = -1;
    }
}

bool TelemLink::init()
{
    // O_NONBLOCK 保证 ::read 时如果不带数据立刻返回，不会阻塞算法主循环
    _fd = ::open(_device_path, O_RDWR | O_NOCTTY | O_NONBLOCK);

    if (_fd < 0) {
        PX4_ERR("TelemLink: open %s failed: %d", _device_path, errno);
        return false;
    }

    struct termios uart_config {};
    if (tcgetattr(_fd, &uart_config) != 0) {
        PX4_ERR("TelemLink: tcgetattr failed");
        ::close(_fd);
        _fd = -1;
        return false;
    }

    cfmakeraw(&uart_config);

    // 设置波特率 115200 8N1
    cfsetispeed(&uart_config, B115200);
    cfsetospeed(&uart_config, B115200);

    uart_config.c_cflag &= ~PARENB;   // no parity
    uart_config.c_cflag &= ~CSTOPB;   // 1 stop bit
    uart_config.c_cflag &= ~CSIZE;
    uart_config.c_cflag |= CS8;       // 8 data bits
    uart_config.c_cflag |= (CLOCAL | CREAD);

#ifdef CRTSCTS
    uart_config.c_cflag &= ~CRTSCTS;  // no flow control
#endif
    uart_config.c_iflag &= ~(IXON | IXOFF | IXANY);

    if (tcsetattr(_fd, TCSANOW, &uart_config) != 0) {
        PX4_ERR("TelemLink: tcsetattr failed");
        ::close(_fd);
        _fd = -1;
        return false;
    }

    PX4_INFO("TelemLink: initialized on %s", _device_path);
    return true;
}

void TelemLink::send_string(const char *str)
{
    if (_fd >= 0 && str != nullptr) {
        ssize_t bytes_written = ::write(_fd, str, strlen(str));
        if (bytes_written < 0) {
            // 如果串口断开了或者写入失败，什么都不做，直接静默忽略
            // （满足了编译器“检查返回值”的强迫症）
        }
    }
}

void TelemLink::update()
{
    if (_fd < 0) {
        return;
    }

    uint8_t ch = 0;
    int n = ::read(_fd, &ch, 1);

    while (n > 0) {
        if (ch == '\r') {
            // 忽略回车符
        } else if (ch == '\n') {
            // 遇到换行符，认为一行结束
            _rx_line[_rx_len] = '\0';
            if (_rx_len > 0) {
                handle_uart_line(_rx_line);
            }
            _rx_len = 0; // 重置游标
        } else {
            // 存入缓冲区
            if (_rx_len < (int)sizeof(_rx_line) - 1) {
                _rx_line[_rx_len++] = (char)ch;
            } else {
                // 缓冲区满（异常过长的残缺帧），直接丢弃，防止内存越界
                _rx_len = 0;
            }
        }

        // 继续读下一个字节，直到缓冲区空（返回 -1）
        n = ::read(_fd, &ch, 1);
    }
}

// ... 前面的 init 和 update 代码保持不变 ...

void TelemLink::handle_uart_line(const char *line)
{
    if (line == nullptr || line[0] == '\0') {
        return;
    }

    // 解析你真正需要的数据
    float rpm = 0.0f;
    if (sscanf(line, "RPM:%f", &rpm) == 1) {

        // 提取成功！存入内部结构体并打上“已更新”标签
        _latest_data.rpm = rpm;
        _latest_data.updated = true;

        // 回复 ACK
        char tx_buf[64];
        snprintf(tx_buf, sizeof(tx_buf), "ACK:RPM:%.2f\r\n", (double)rpm);
        send_string(tx_buf);
    }
    // ... 可以继续增加其他格式的 sscanf 解析 ...
}

// 新增提取接口的实现
bool TelemLink::get_latest_data(TelemData &out_data)
{
    if (_latest_data.updated) {
        out_data = _latest_data;    // 拷贝给外部
        _latest_data.updated = false; // 阅后即焚（清除标志位）
        return true;
    }
    return false;
}
