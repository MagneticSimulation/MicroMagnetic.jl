using Luxor
using Colors

function draw_micromagnetic_logo(output_path="logo.svg")
    julia_colors = [ 
        Luxor.julia_purple, 
        Luxor.julia_green, 
        Luxor.julia_red, 
        Luxor.julia_blue
    ]
    
    width, height = 500, 500
    Drawing(width, height, output_path)
    #background("white")
    origin()
    
    square_size = 200
    rotate(π/4)
    
    corners = [
        Point(-square_size/2, -square_size/2),
        Point(square_size/2, -square_size/2),
        Point(square_size/2, square_size/2),
        Point(-square_size/2, square_size/2)
    ]
    
    circle_radius = 80
    arrow_length = 80
    
    arrow_directions = [-π/4, π/4, 3π/4, 5π/4]
    
    for (i, corner) in enumerate(corners)
        setcolor(julia_colors[i])
        circle(corner, circle_radius, :fill)
        setcolor("black")
        angle = arrow_directions[i]
        start_point = corner - 0.5*arrow_length * Point(cos(angle), sin(angle))
        end_point = start_point + arrow_length * Point(cos(angle), sin(angle))
        arrow(start_point, end_point, linewidth=16, arrowheadlength=36, arrowheadangle=pi/6)
    end
    
    finish()
    println("Logo saved to: ", output_path)
end

draw_micromagnetic_logo("logo.svg")