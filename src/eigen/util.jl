struct Epsilon <: AbstractFloat
    id::Int
    value::Float64
end

struct AddExpr <: AbstractFloat
    a::AbstractFloat
    b::AbstractFloat
end

# Promotion rules, we don't have to define them
Base.promote_rule(::Type{Epsilon}, ::Type{Epsilon}) = Epsilon
Base.promote_rule(::Type{Epsilon}, ::Type{Real}) = Epsilon
Base.promote_rule(::Type{Real}, ::Type{Epsilon}) = Epsilon
Base.promote_rule(::Type{AddExpr}, ::Type{Real}) = AddExpr
Base.promote_rule(::Type{Real}, ::Type{AddExpr}) = AddExpr
Base.promote_rule(::Type{AddExpr}, ::Type{Epsilon}) = AddExpr
Base.promote_rule(::Type{Epsilon}, ::Type{AddExpr}) = AddExpr

# multiplication 
Base.:*(a::Epsilon, b::Epsilon) = 0
Base.:*(a::Real, b::Epsilon) = a == 0 ? 0 : Epsilon(b.id, a * b.value)
Base.:*(a::Epsilon, b::Real) = b * a
Base.:*(a::Real, b::AddExpr) = a == 0 ? 0 : AddExpr(a * b.a, a * b.b)
Base.:*(a::AddExpr, b::Real) = b * a
Base.:*(a::Epsilon, b::AddExpr) = AddExpr(a * b.a, a * b.b)
Base.:*(a::AddExpr, b::AddExpr) = AddExpr(a.a * b, a.b * b)
Base.:*(a::AddExpr, b::Epsilon) = b * a

# add
function Base.:+(a::Epsilon, b::Epsilon)
    return a.id == b.id ? Epsilon(a.id, a.value + b.value) : AddExpr(a, b)
end
Base.:+(a::Real, b::Epsilon) = a == 0 ? b : AddExpr(a, b)
Base.:+(a::Epsilon, b::Real) = b + a
Base.:+(a::Real, b::AddExpr) = a == 0 ? b : AddExpr(a, b)
Base.:+(a::AddExpr, b::Real) = b + a
Base.:+(expr::AddExpr, eps::Epsilon) = AddExpr(expr, eps)
Base.:+(eps::Epsilon, expr::AddExpr) = expr + eps
Base.:+(e1::AddExpr, e2::AddExpr) = AddExpr(e1, e2)

Base.:-(b::Epsilon) = Epsilon(b.id, -b.value)
Base.:-(b::AddExpr) = AddExpr(-b.a, -b.b)

# subtraction
function Base.:-(a::Epsilon, b::Epsilon)
    return a.id == b.id ? Epsilon(a.id, a.value - b.value) : AddExpr(a, -b)
end
Base.:-(a::Real, b::Epsilon) = a + (-b)
Base.:-(a::Epsilon, b::Real) = a + (-b)
Base.:-(a::Real, b::AddExpr) = a + (-b)
Base.:-(a::AddExpr, b::Real) = AddExpr(a.a - b, a.b)
Base.:-(a::Epsilon, b::AddExpr) = AddExpr(a, -b)
Base.:-(a::AddExpr, b::Epsilon) = AddExpr(a.a - b, a.b)
Base.:-(e1::AddExpr, e2::AddExpr) = AddExpr(e1, -e2)

# display 
function Base.show(io::IO, expr::AbstractFloat)
    if expr isa Epsilon
        if expr.value >= 0
            print(io, expr.value, "ε_", expr.id)
        else
            print(io, "(", expr.value, "ε_", expr.id, ")")
        end
    elseif expr isa AddExpr
        print(io, "(", expr.a, " + ", expr.b, ")")
    else
        print(io, expr)
    end
end

function collect_terms(expr::AbstractFloat)
    terms = Dict{Int,AbstractFloat}()
    function collect!(expr)
        if expr isa Epsilon
            terms[expr.id] = get(terms, expr.id, 0.0) + expr.value
        elseif expr isa AddExpr
            collect!(expr.a)
            collect!(expr.b)
        elseif expr isa AbstractFloat
            terms[0] = get(terms, 0, 0.0) + expr
        end
    end
    collect!(expr)
    return terms
end


function simplify(expr::AbstractFloat)
    terms = collect_terms(expr)

    result = haskey(terms, 0) ? terms[0] : 0

    sorted_keys = sort(collect(keys(terms)))
    for key in sorted_keys
        if key != 0
            result = result + Epsilon(key, terms[key])
        end
    end

    return result
end

raw"""
R z = m0
where R =
 \begin{bmatrix}
 \cos \phi \cos \theta & -\sin \phi & \sin \theta \cos \phi \\
 \sin \phi \cos \theta &  \cos \phi & \sin \theta \sin \phi \\
 -\sin \theta & 0 &  \cos \theta \\
 \end{bmatrix}

and
    m_x = \sin \theta \cos \phi
    m_y = \sin \theta \sin \phi
    m_z = \cos \theta
"""
function rotation_matrix(mx, my, mz)
    norm = sqrt(mx^2 + my^2 + mz^2)
    A = zeros(3,3)
    if norm < eps()
        A[1,1] = 1
        A[2,2] = 1
        A[3,3] = 1
    else
        theta = acos(mz / norm)
        phi = atan(my, mx)
        A[1, 1]  = cos(phi) * cos(theta)
        A[1, 2] = -sin(phi) 
        A[1, 3] = sin(theta) * cos(phi)
        A[2, 1] = sin(phi) * cos(theta)
        A[2, 2] = cos(phi)
        A[2, 3] = sin(theta) * sin(phi)
        A[3, 1] = -sin(theta)
        A[3, 2] = 0
        A[3, 3] = cos(theta)
    end
    return A
end

"""
R^-1 = R^T where R is given by rotation_matrix

R^-1 m0 = z
"""
function rotation_matrix_inverse(mx, my, mz)
    A = rotation_matrix(mx, my, mz)
    return transpose(A)
end


@inline function cross_x(x1, x2, x3, y1, y2, y3)
    return -x3*y2 + x2*y3
end

@inline function cross_y(x1, x2, x3, y1, y2, y3)
    return x3*y1 - x1*y3
end

@inline function cross_z(x1, x2, x3, y1, y2, y3)
    return -x2*y1 + x1*y2
end