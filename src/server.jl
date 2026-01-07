using HTTP
using JSON
using Dates

export start_server, serve, gui

const DEFAULT_HOST = "0.0.0.0"
const DEFAULT_PORT = 10056
global current_port = DEFAULT_PORT
# global variables for managing single connection
global_sim = nothing
global_client = nothing
# Message handlers
message_handlers = Dict{String, Function}()

function run_code(code::String)
    stdout_pipe = Pipe()
    stderr_pipe = Pipe()
    
    success = true
    result = nothing
    try
        result = redirect_stdio(stdout=stdout_pipe, stderr=stderr_pipe) do
            eval(Meta.parseall(code))
        end
    catch e
        success = false
        print(stderr_pipe.in, sprint(showerror, e))
    finally
        close(stdout_pipe.in)
        close(stderr_pipe.in)
    end
    
    stdout_str = read(stdout_pipe.out, String)
    stderr_str = read(stderr_pipe.out, String)
    
    close(stdout_pipe.out)
    close(stderr_pipe.out)
    
    if success && !isempty(stderr_str)
        stdout_str = stdout_str * (isempty(stdout_str) ? "" : "\n") * stderr_str
        stderr_str = ""
    end
        
    return result, success, stdout_str, stderr_str
end

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


"""
Log a message with timestamp
"""
function log_message(message::String)
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
    @info "[$timestamp] $message"
end

function register_handler(handler::Function, message_type::String)
    message_handlers[message_type] = handler
end

"""
Initialize default message handlers
"""
function init_default_handlers()
    
    register_handler("run_code") do ws, data
        code = get(data, "code", "")
        command_id = get(data, "id", "")
        description = get(data, "desc", "")
        
        log_message("Received run_code command (ID: $command_id, Desc: $description)")
        
        result, success, stdout_str, stderr_str = run_code(code)
        response = Dict(
            "id" => command_id,
            "success" => success,
            "stdout" => stdout_str,       
            "stderr" => stderr_str,
            "type" => description
        )

        if isa(result, FDMesh)
            response["type"] = "fd_mesh_data"
            response["fd_mesh_data"] = Dict(
                "nx" => result.nx,
                "ny" => result.ny,
                "nz" => result.nz,
                "dx" => result.dx*1e9,
                "dy" => result.dy*1e9,
                "dz" => result.dz*1e9,
            )
        elseif isa(result, MicroSim)
           global global_sim = result
        end

        if description == "ms" && global_sim !== nothing
            response["type"] = "Ms_data"
            response["Ms_data"] = global_sim.mu0_Ms
        elseif (description == "m0" || description == "relax") && global_sim !== nothing
            response["type"] = "m_data"
            response["m_data"] = Array(global_sim.spin)
        end

        send_message(ws, "run_code_response", response)
    end
    
    register_handler("stop_server") do ws, client_id, data
        log_message("Received stop_server command from client $client_id")
        
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
    # Since we only allow one connection, simplify broadcasting
    if global_client !== nothing
        send_message(global_client, message_type, data)
    end
end

"""
Handle WebSocket message
"""
function handle_message(ws, message::String)
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
    # Check if there's already an active connection
    if global_client !== nothing
        log_message("Connection rejected: server already has an active client connection")
        # Send error message and close connection
        try
            send_message(ws, "error", Dict(
                "message" => "Server already has an active client connection. Only one connection is allowed at a time."
            ))
            close(ws)
        catch e
            @warn "Error sending rejection message: $e"
        end
        return
    end
    
    global global_client = ws
    log_message("WebSocket client connected: $client_id")
    
    try     
        # Listen for messages
        for msg in ws
            handle_message(ws, msg)
        end
        
    catch e
        if e isa HTTP.WebSockets.WebSocketError
            log_message("WebSocket client disconnected: $client_id")
        else
            log_message("WebSocket error: $e")
        end
    finally
        global global_client = nothing
        global global_sim = nothing
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
            # active_connections[client_id] = ws
            handle_websocket(ws, client_id)
        end
        
    elseif req.method == "GET"
        # HTTP GET request
        handle_get_request(http, req, message_handlers)
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
function handle_get_request(http, req, message_handlers)
    if req.target == "/status"
        # Return JSON status
        status = Dict(
            "status" => "running",
            "port" => current_port,
            "active_connections" => global_client !== nothing ? 1 : 0,
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
        
        full_path = joinpath(@__DIR__, "../gui", file_path)
        
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
                # Client closed connection
            end
        else
            # File not found
            HTTP.setstatus(http, 404)
            HTTP.setheader(http, "Content-Type" => "text/plain")
            HTTP.startwrite(http)
            HTTP.write(http, "File not found")
        end
    end
end