using HTTP
using JSON
using Dates

export start_server, serve, gui

const DEFAULT_HOST = "0.0.0.0"
const DEFAULT_PORT = 10056
global current_port = DEFAULT_PORT

"""
Decode URI component
"""
function decodeURIComponent(uri::AbstractString)
    # Replace %20 with space
    uri = replace(uri, "%20" => " ")
    # Replace other %xx sequences
    uri = replace(uri, r"%([0-9a-fA-F]{2})" => s -> string(Char(parse(UInt8, s, base=16))))
    return uri
end

# Global variables
global current_port = 0
global active_connections = Dict{String, HTTP.WebSockets.WebSocket}()
global message_handlers = Dict{String, Function}()

"""
Log a message with timestamp
"""
function log_message(message::String)
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
    @info "[$timestamp] $message"
end

"""
Message handler type
"""
struct MessageHandler
    message_type::String
    handler::Function
end

function register_handler(handler::Function, message_type::String)
    message_handlers[message_type] = handler
end

"""
Initialize default message handlers
"""
function init_default_handlers()
    # Register default handlers
    register_handler("ping") do ws, data
        send_message(ws, "pong", Dict("timestamp" => now()))
    end
    
    register_handler("echo") do ws, data
        send_message(ws, "echo_response", Dict("message" => get(data, "message", "")))
    end
    
    register_handler("get_status") do ws, data
        status = Dict( 
            "server_time" => now(), 
            "connections" => length(active_connections), 
            "status" => "running" 
        ) 
        send_message(ws, "status", status) 
    end
    
    register_handler("simulation_start") do ws, data
        # Integration with Micromagnetic.jl simulation logic can be added here
        println("Starting simulation: ", data)
        send_message(ws, "simulation_started", Dict(
            "simulation_id" => randstring(10),
            "start_time" => now()
        ))
    end
    
    register_handler("simulation_stop") do ws, data
        println("Stopping simulation: ", data)
        send_message(ws, "simulation_stopped", Dict("stop_time" => now()))
    end
    
    register_handler("get_data") do ws, data
        # Sample data return - replace with Micromagnetic.jl data in actual use
        sample_data = Dict(
            "time_points" => collect(0:0.1:10),
            "magnetization" => rand(101, 3),
            "energy" => rand(101),
            "parameters" => data
        )
        send_message(ws, "data_response", sample_data)
    end
    
    # Maintain backward compatibility
    register_handler("command") do ws, data
        command = get(data, "command", "")
        command_id = get(data, "id", "")
        
        log_message("Received command: $command (ID: $command_id)")
        response = Dict(
            "id" => command_id,
            "command" => command,
            "success" => true,
            "result" => Dict("status" => "processed", "server_time" => now())
        )
        send_message(ws, "command_response", response)
    end
    
    register_handler("heartbeat") do ws, data
        send_message(ws, "heartbeat_response", Dict("timestamp" => now()))
    end
    
    register_handler("stop_server") do ws, data
        # Get client_id from the active_connections (find the key for this ws)
        client_id = "unknown"
        for (id, connection) in active_connections
            if connection === ws
                client_id = id
                break
            end
        end
        
        log_message("Received stop_server command from client $client_id")
        
        # Instead of stopping the entire server, just close the client's connection
        # This way, other clients are not affected
        
        try
            # Send confirmation message to the client
            send_message(ws, "server_stopping", Dict(
                "message" => "Your connection is being closed",
                "timestamp" => now()
            ))
            
            # Close the WebSocket connection
            close(ws)
            
            log_message("Closed connection for client $client_id")
        catch e
            log_message("Error closing connection for client $client_id: $e")
        end
    end
end

"""
Send message to client
"""
function send_message(ws, message_type::String, data::Dict=Dict())
    message = Dict(
        "type" => message_type,
        "data" => data,
        "timestamp" => string(now())
    )
    try
        HTTP.WebSockets.send(ws, JSON.json(message))
        return true
    catch e
        @error "Failed to send message" exception=e
        return false
    end
end

"""
Broadcast message to all connected clients
"""
function broadcast(message_type::String, data::Dict=Dict())
    for (client_id, ws) in active_connections
        send_message(ws, message_type, data)
    end
end

"""
Handle WebSocket message
"""
function handle_message(ws, client_id::String, message::String)
    try
        # Parse JSON message
        msg = JSON.parse(message)
        
        # Check message format
        if !haskey(msg, :type)
            error_msg = Dict("error" => "Message missing 'type' field")
            send_message(ws, "error", error_msg)
            return
        end
        
        message_type = string(msg.type)
        data = haskey(msg, :data) ? msg.data : Dict()
        
        # Find and execute handler
        if haskey(message_handlers, message_type)
            try
                message_handlers[message_type](ws, data)
            catch e
                @error "Failed to handle message" message_type exception=e
                send_message(ws, "error", Dict(
                    "message" => "Failed to handle message: $message_type",
                    "error" => string(e)
                ))
            end
        else
            @warn "Unknown message type: $message_type"
            send_message(ws, "error", Dict(
                "message" => "Unknown message type: $message_type"
            ))
        end
        
    catch e
        @error "Failed to parse message" message exception=e
        send_message(ws, "error", Dict(
            "message" => "Message format error",
            "error" => string(e)
        ))
    end
end

"""
WebSocket connection handler
"""
function handle_websocket(ws, client_id::String)
    log_message("WebSocket client connected: $client_id")
    
    try
        # Send connection confirmation
        send_message(ws, "connected", Dict(
            "client_id" => client_id,
            "server_info" => "Micromagnetic.jl WebSocket Server",
            "timestamp" => now()
        ))
        
        # Listen for messages
        for msg in ws
            handle_message(ws, client_id, msg)
        end
        
    catch e
        if e isa HTTP.WebSockets.WebSocketError
            log_message("WebSocket client disconnected: $client_id")
        else
            log_message("WebSocket error: $e")
        end
    finally
        # Clean up connection
        delete!(active_connections, client_id)
        log_message("Connection closed: $client_id")
    end
end

"""
Generate client ID
"""
function generate_client_id()
    return "client_" * randstring(8)
end

"""
HTTP request handler
"""
function handle_http_request(http::HTTP.Stream)
    req = http.message
    
    if HTTP.WebSockets.isupgrade(req)
        # WebSocket upgrade request
        # Try to get session ID from headers or query parameters
        session_id = nothing
        
        # Check headers for session ID
        # Find session ID in headers vector
        session_id = nothing
        for (key, value) in req.headers
            if key == "X-Session-ID"
                session_id = value
                break
            end
        end
        
        # Check query parameters for session ID
        if session_id === nothing && !isempty(req.target)
            query_start = findfirst('?', req.target)
            if query_start !== nothing
                query_str = req.target[query_start+1:end]
                query_params = Dict{String, String}()
                for param in split(query_str, '&')
                    if contains(param, '=')
                        key, value = split(param, '=', limit=2)
                        query_params[decodeURIComponent(key)] = decodeURIComponent(value)
                    end
                end
                if haskey(query_params, "session_id")
                    session_id = query_params["session_id"]
                end
            end
        end
        
        # Use session ID as client ID, or generate a new one if not provided
        client_id = if session_id !== nothing
            session_id
        else
            generate_client_id()
        end
        
        HTTP.WebSockets.upgrade(http) do ws
            active_connections[client_id] = ws
            handle_websocket(ws, client_id)
        end
        
    elseif req.method == "GET"
        # HTTP GET request
        handle_get_request(http, req, active_connections, message_handlers)
    else
        # Method not allowed
        HTTP.setstatus(http, 405)
        HTTP.setheader(http, "Content-Type" => "text/plain")
        HTTP.startwrite(http)
        HTTP.write(http, "Method not allowed")
    end
end

"""
Start MicroMagnetic GUI server
"""
function start_server(port::Int=DEFAULT_PORT, host::String=DEFAULT_HOST)
    max_attempts = 5
    
    for attempt in 1:max_attempts
        current_attempt_port = port + (attempt - 1)
        
        try
            log_message("Starting MicroMagnetic server on $host:$current_attempt_port (attempt $attempt/$max_attempts)")
            println("Server address: http://$host:$current_attempt_port")
            println("Press Ctrl+C to stop server")
            
            # Initialize default message handlers
            init_default_handlers()
            
            # Update global current_port variable
            global current_port = current_attempt_port
            
            # Create a server object that can be stopped
            server = nothing
            
            # Try to start the server
            server = HTTP.listen(host, current_attempt_port) do http::HTTP.Stream
                handle_http_request(http)
            end
            
            # If we get here, the server started successfully and will run until interrupted
            return
            
        catch e
            if e isa InterruptException
                log_message("Server stopped by user")
                return
            elseif occursin("address already in use", lowercase(string(e))) || 
                   occursin("port is already in use", lowercase(string(e))) ||
                   occursin("listen: address already in use", lowercase(string(e)))
                log_message("Port $current_attempt_port is already in use, trying port $(current_attempt_port + 1)...")
                # Continue to next iteration to try the next port
            else
                log_message("Server error: $e")
                # For other errors, stop trying
                rethrow(e)
            end
        end
    end
    
    # If we've tried all attempts and failed
    error("Failed to start server after $max_attempts attempts. All ports from $port to $(port + max_attempts - 1) are in use.")
end

"""
Start MicroMagnetic GUI (convenience function)
"""
function gui(port::Int=DEFAULT_PORT, host::String=DEFAULT_HOST)
    println("Starting MicroMagnetic GUI...")
    println("GUI will be available at: http://$host:$port")
    println("Press Ctrl+C to exit")
    start_server(port, host)
end

"""
Handle HTTP GET requests
"""
function handle_get_request(http, req, active_connections, message_handlers)
    if req.target == "/status"
        # Return JSON status
        status = Dict(
            "status" => "running",
            "port" => current_port,
            "active_connections" => length(active_connections),
            "server_time" => string(now()),
            "message_handlers" => collect(keys(message_handlers))
        )
        
        try
            HTTP.setstatus(http, 200)
            HTTP.setheader(http, "Content-Type" => "application/json")
            HTTP.startwrite(http)
            HTTP.write(http, JSON.json(status))
        catch e
            # Client closed connection
        end
    else
        # Serve static files (including root path "/")
        file_path = lstrip(req.target, '/')
        if file_path == ""
            file_path = "index.html"
        end
        
        full_path = joinpath(dirname(@__FILE__), "js", file_path)
        
        if isfile(full_path)
            # Determine MIME type
            ext = lowercase(splitext(file_path)[2])
            mime_type = if ext == ".html"
                "text/html"
            elseif ext == ".js"
                "application/javascript"
            elseif ext == ".css"
                "text/css"
            elseif ext == ".json"
                "application/json"
            else
                "application/octet-stream"
            end
            
            try
                content = read(full_path)
                HTTP.setstatus(http, 200)
                HTTP.setheader(http, "Content-Type" => mime_type)
                HTTP.startwrite(http)
                HTTP.write(http, content)
            catch e
                
            end
        else
            HTTP.setstatus(http, 404)
            HTTP.setheader(http, "Content-Type" => "text/plain")
            HTTP.startwrite(http)
            HTTP.write(http, "File not found")
        end
    end
end