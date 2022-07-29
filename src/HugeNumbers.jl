module HugeNumbers

export Huge, HugeFloat16, HugeFloat32, HugeFloat64, hugen, invhugen

issmall(x) = abs(x) < one(x)

"""
    hugen(x)

Compute the `hugen` function at `x`.

It is defined as:
- `exp(x - 1)` if `x ≥ 1`
- `exp(1 - 1/x)` if `x ≥ 0`
- `-exp(1 + 1/x)` if `x ≥ -1`
- `-exp(-x - 1)` otherwise

It has these nice arithmetic properties:
- It is an increasing, continuous, differentiable, invertible function on the real numbers.
- It grows exponentially for large `x` and shrinks exponentially for small `x`.
- `hugen(x) = x` for `x` in `-Inf`, `-1`, `0`, `1`, `Inf`.
- `-hugen(x) = hugen(-x)`
- `1/hugen(x) = hugen(1/x)`
- Its derivative is `1` at `±1`, and `0` at `0`.

The inverse function is [`invhugen`](@ref).
"""
function hugen(x)
    if signbit(x)
        if issmall(x)
            return -exp(one(x) + inv(x))
        else
            return -exp(-x - one(x))
        end
    else
        if issmall(x)
            return exp(one(x) - inv(x))
        else
            return exp(x - one(x))
        end
    end
end

"""
    invhugen(x)

The inverse function of [`hugen`](@ref).

That is, `invhugen(hugen(x)) == x` and `hugen(invhugen(x)) == x` (up to floating point
precision).
"""
function invhugen(x)
    if signbit(x)
        y = log(-x)
        if issmall(x)
            return inv(y - one(y))
        else
            return -y - one(y)
        end
    else
        y = log(x)
        if issmall(x)
            return inv(one(y) - y)
        else
            return y + one(y)
        end
    end
end

### CONSTRUCTION / CONVERSION / PROMOTION

"""
    Huge{[T]}(x)

Convert `x` to a `Huge` number, which is represented by `invhugen(x)`.

If you already know `ix = invhugen(x)`, then you may call `hugen(Huge, ix)` instead.
"""
struct Huge{T<:Real} <: Real
    inv::T
    global hugen(::Type{Huge{T}}, x::T) where {T<:Real} = new{T}(x)
end

const HugeFloat16 = Huge{Float16}
const HugeFloat32 = Huge{Float32}
const HugeFloat64 = Huge{Float64}

"""
    hugen(Huge, x)

Compute `hugen(x)` but the result is stored exactly as a `Huge`.
"""
hugen(::Type{Huge{T}}, x::Real) where {T<:Real} = hugen(Huge{T}, convert(T, x))
hugen(::Type{Huge}, x::T) where {T<:Real} = hugen(Huge{T}, x)

invhugen(x::Huge) = x.inv
hugen(x::Huge) = hugen(Huge, hugen(invhugen(x)))

Base.convert(::Type{Huge{T}}, x::Huge{T}) where {T} = x
Base.convert(::Type{Huge{T}}, x::Real) where {T} = hugen(Huge{T}, invhugen(x))

Base.convert(::Type{Huge}, x::Huge) = x
Base.convert(::Type{Huge}, x::Real) = hugen(Huge, invhugen(x))

Huge(x::Real) = convert(Huge, x)

Huge{T}(x::Real) where {T} = convert(Huge{T}, x)

Base.AbstractFloat(x::Huge) = AbstractFloat(hugen(AbstractFloat(invhugen(x))))

Base.BigFloat(x::Huge) = BigFloat(hugen(BigFloat(invhugen(x))))

Base.promote_rule(::Type{Huge{T}}, ::Type{Huge{S}}) where {T<:Real, S<:Real} = Huge{promote_type(T, S)}
Base.promote_rule(::Type{Huge{T}}, ::Type{S}) where {T<:Real, S<:Real} = promote_type(Huge{T}, Huge{S})


### COMPARISONS / PREDICATES

Base.signbit(x::Huge) = signbit(invhugen(x))

Base.sign(x::Huge) = hugen(Huge, sign(invhugen(x)))

Base.abs(x::Huge) = hugen(Huge, abs(invhugen(x)))

Base.isfinite(x::Huge) = isfinite(invhugen(x))

Base.isinf(x::Huge) = isinf(invhugen(x))

Base.isnan(x::Huge) = isnan(invhugen(x))

Base.iszero(x::Huge) = iszero(invhugen(x))

Base.isone(x::Huge) = isone(invhugen(x))

Base.isless(x::Huge, y::Huge) = isless(invhugen(x), invhugen(y))

Base.isequal(x::Huge, y::Huge) = isequal(invhugen(x), invhugen(y))

Base.:(==)(x::Huge, y::Huge) = isless(invhugen(x), invhugen(y))

Base.:(<)(x::Huge, y::Huge) = invhugen(x) < invhugen(y)

Base.cmp(x::Huge, y::Huge) = cmp(invhugen(x), invhugen(y))

Base.zero(::Type{Huge}) = zero(Huge{Int})
Base.zero(::Type{Huge{T}}) where {T<:Real} = hugen(Huge, zero(T))

Base.one(::Type{Huge}) = one(Huge{Int})
Base.one(::Type{Huge{T}}) where {T<:Real} = hugen(Huge, one(T))

# We hash the inner value, which for free means that hash(Huge(x)) == hash(x) for x in
# [-Inf, -1, 0, 1, Inf], which are also the only values likely to overlap with other types.
Base.hash(x::Huge, h::UInt) = hash(invhugen(x), h)

Base.typemin(::Type{Huge{T}}) where {T} = hugen(Huge, typemin(T))

Base.typemax(::Type{Huge{T}}) where {T} = hugen(Huge, typemax(T))

Base.nextfloat(x::Huge) = hugen(Huge, nextfloat(invhugen(x)))
Base.nextfloat(x::Huge, n::Integer) = hugen(Huge, nextfloat(invhugen(x), n))

Base.prevfloat(x::Huge) = hugen(Huge, prevfloat(invhugen(x)))
Base.prevfloat(x::Huge, n::Integer) = hugen(Huge, prevfloat(invhugen(x), n))

### TYPES

Base.big(::Type{Huge}) = Huge{BigFloat}
Base.big(::Type{Huge{T}}) where {T} = Huge{big(T)}
Base.big(x::Huge) = convert(big(typeof(x)), x)

Base.widen(::Type{Huge{T}}) where {T} = Huge{widen(T)}
Base.widen(x::Huge) = convert(widen(typeof(x)), x)

Base.float(::Type{Huge{T}}) where {T} = Huge{float(T)}

### ARITHMETIC

Base.:(-)(x::Huge) = hugen(Huge, -invhugen(x))

Base.inv(x::Huge) = hugen(Huge, inv(invhugen(x)))

function Base.:(*)(x::Huge{T}, y::Huge{T}) where {T}
    ix = invhugen(x)
    iy = invhugen(y)
    if signbit(ix)
        if signbit(iy)
            ir = _mul(-ix, -iy)
        else
            ir = -_mul(-ix, iy)
        end
    else
        if signbit(iy)
            ir = -_mul(ix, -iy)
        else
            ir = _mul(ix, iy)
        end
    end
    return hugen(Huge, ir)
end

Base.:(/)(x::Huge{T}, y::Huge{T}) where {T} = x * inv(y)

function _mul(ix, iy)
    if issmall(ix)
        if issmall(iy)
            inv(__mul(inv(ix), inv(iy)))
        else
            __div(iy, inv(ix))
        end
    else
        if issmall(iy)
            __div(ix, inv(iy))
        else
            __mul(ix, iy)
        end
    end
end

function __mul(ix, iy)
    iz = ix + iy
    return iz - one(iz)
end

function __div(ix, iy)
    iz = ix - iy
    if signbit(iz)
        return inv(-iz + one(iz))
    else
        return iz + one(iz)
    end
end

function Base.:(+)(x::Huge{T}, y::Huge{T}) where {T}
    ix = invhugen(x)
    iy = invhugen(y)
    if signbit(x)
        if signbit(y)
            ir = -_add(-ix, -iy)
        else
            ir = _sub(iy, -ix)
        end
    else
        if signbit(y)
            ir = _sub(ix, -iy)
        else
            ir = _add(ix, iy)
        end
    end
    return hugen(Huge, ir)
end

Base.:(-)(x::Huge{T}, y::Huge{T}) where {T} = x + (-y)

function _add(ix, iy)
    if issmall(ix)
        if issmall(iy)
            if ix < iy
                t = inv(iy) - log1p(exp(inv(iy) - inv(ix)))
            else
                t = inv(ix) - log1p(exp(inv(ix) - inv(iy)))
            end
            return t < one(t) ? one(t) + one(t) - t : inv(t)
        else
            return iy + log1p(exp((one(ix) - inv(ix)) - (iy - one(iy))))
        end
    else
        if issmall(iy)
            return ix + log1p(exp((one(iy) - inv(iy)) - (ix - one(ix))))
        else
            if ix < iy
                return iy + log1p(exp(ix - iy))
            else
                return ix + log1p(exp(iy - ix))
            end
        end
    end
end

function _sub(ix, iy)
    if ix < iy
        return -__sub(iy, ix)
    else
        return __sub(ix, iy)
    end
end

function __sub(ix, iy)
    if issmall(ix)
        @assert issmall(iy) || isnan(iy)
        return inv(inv(ix) - log1p(-exp(inv(ix) - inv(iy))))
    else
        if issmall(iy)
            ir = ix + log1p(-exp((one(iy) - inv(iy)) - (ix - one(ix))))
        else
            ir = ix + log1p(-exp(iy - ix))
        end
        return ir < one(ir) ? inv(one(ir) + one(ir) - ir) : ir
    end
end

function Base.log(x::Huge)
    ix = invhugen(x)
    if signbit(ix) && !iszero(ix) && !isnan(ix)
        throw(DomainError(x))
    elseif issmall(ix)
        return one(ix) - inv(ix)
    else
        return ix - one(ix)
    end
end

function Base.log2(x::Huge)
    logx = log(x)
    return logx / log(oftype(logx, 2))
end

function Base.log10(x::Huge)
    logx = log(x)
    return logx / log(oftype(logx, 2))
end

function Base.exp(::Type{Huge}, x::Real)
    if signbit(x)
        ir = inv(one(x) - x)
    else
        ir = x + inv(inv(one(x)))
    end
    return hugen(Huge, ir)
end

function Base.exp(::Type{Huge{T}}, x::Real) where {T}
    if signbit(x)
        ir = inv(one(T) - convert(T, x))
    else
        ir = convert(T, x) + one(T)
    end
    return hugen(Huge{T}, ir)
end


## IO

function Base.show(io::IO, x::Huge)
    print(io, "hugen(")
    if get(io, :typeinfo, Any) != typeof(x)
        show(io, typeof(x))
        print(io, ", ")
    end
    show(io, invhugen(x))
    print(io, ")")
end

function Base.write(io::IO, x::Huge)
    write(io, invhugen(x))
end

function Base.read(io::IO, ::Type{Huge{T}}) where {T}
    hugen(Huge{T}, read(io, T))
end

end # module
