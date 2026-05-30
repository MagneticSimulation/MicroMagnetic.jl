# GUI

MicroMagnetic provides a web-based graphical user interface for visualizing simulation results, executing code, and controlling simulations.

## Demo Video

Check out this demo video: [MicroMagnetic GUI Demo](https://meeting.tencent.com/crm/NQGpVpWW7c)

## Starting the GUI

```julia
using MicroMagnetic

# Start GUI (default port 10056)
gui()

# Start GUI accessible from LAN
gui(lan=true)

# Start GUI on custom port
gui(port=8080)
```

Or use the `start_server` function:

```julia
using MicroMagnetic

# Synchronous mode
start_server(async=false)

# Asynchronous mode (default)
start_server(port=10056, host=nothing, async=true, lan=true)
```

## Important Notes

1. **Shared Global State**: Code executed in the Julia terminal and code executed in the GUI share the same global state. This means variables created in either environment are visible to both.

2. **One-to-One Mapping**: Each terminal can only connect to one GUI web page. The GUI communicates with the Julia session running in the terminal that started it.

3. **Multiple GUIs**: To run multiple GUI instances:
   - Start multiple Julia terminals
   - Each terminal runs on a different port
   - Connect each browser to its corresponding port

Example:
```julia
# Terminal 1
gui(port=10056)  # Browser: http://localhost:10056

# Terminal 2  
gui(port=10057)  # Browser: http://localhost:10057

# Terminal 3
gui(port=10058)  # Browser: http://localhost:10058
```

