# 图形界面

MicroMagnetic 提供了一个基于网页的图形界面，用于可视化模拟结果、执行代码以及控制模拟。

## 演示视频

观看演示视频：[MicroMagnetic GUI 演示](https://meeting.tencent.com/crm/NQGpVpWW7c)

## 启动 GUI

```julia
using MicroMagnetic

# 启动 GUI（默认端口 10056）
gui()

# 启动可供局域网访问的 GUI
gui(lan=true)

# 在自定义端口上启动 GUI
gui(port=8080)
```

也可以使用 `start_server` 函数：

```julia
using MicroMagnetic

# 同步模式
start_server(async=false)

# 异步模式（默认）
start_server(port=10056, host=nothing, async=true, lan=true)
```

## 注意事项

1. **共享全局状态**：在 Julia 终端中执行的代码和在 GUI 中执行的代码共享同一个全局状态。也就是说，在任一环境中创建的变量在两个环境中都可见。

2. **一对一映射**：每个终端只能连接一个 GUI 网页。GUI 与启动它的终端中运行的 Julia 会话通信。

3. **多个 GUI 实例**：要同时运行多个 GUI 实例：
   - 启动多个 Julia 终端
   - 每个终端使用不同的端口
   - 将每个浏览器连接到对应的端口

示例：
```julia
# 终端 1
gui(port=10056)  # 浏览器访问: http://localhost:10056

# 终端 2  
gui(port=10057)  # 浏览器访问: http://localhost:10057

# 终端 3
gui(port=10058)  # 浏览器访问: http://localhost:10058
```
