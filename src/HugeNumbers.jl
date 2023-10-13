module HugeNumbers

export Huge, HugeFloat16, HugeFloat32, HugeFloat64, hexp, hlog

issmall(x) = abs(x) < one(x)

function h(x)
    ax = abs(x)
    ix = inv(ax)
    r1 = ax - one(ax)
    r2 = one(ix) - ix
    r1, r2 = promote(r1, r2)
    return ifelse(ax < one(ax), r2, r1)
end

function invh(x)
    r1 = x + one(x)
    r2 = inv(one(x) - x)
    r1, r2 = promote(r1, r2)
    return ifelse(signbit(x), r2, r1)
end

"""
    hexp(x)

Compute the `hexp` function at `x`.

It is defined as:
- `exp(x - 1)` if `x ≥ 1`
- `exp(1 - 1/x)` if `x ≥ 0`
- `-exp(1 + 1/x)` if `x ≥ -1`
- `-exp(-x - 1)` otherwise
    
It has these nice arithmetic properties:
- It is an increasing, continuous, differentiable, invertible function on the real numbers.
- It grows exponentially for large `x` and shrinks exponentially for small `x`.
- `hexp(x) = x` for `x` in `-Inf`, `-1`, `0`, `1`, `Inf`.
- `-hexp(x) = hexp(-x)`
- `1/hexp(x) = hexp(1/x)`
- Its derivative is `1` at `±1`, and `0` at `0`.

The inverse function is [`hlog`](@ref).
"""
hexp(x) = copysign(exp(h(x)), x)

"""
    hlog(x)

The inverse function of [`hexp`](@ref).
"""
hlog(x) = copysign(invh(log(abs(x))), x)


### CONSTRUCTION / CONVERSION / PROMOTION

"""
    Huge{[T]}(x)

Convert `x` to a `Huge` number, which is represented by `hlog(x)`.

If you already know `ix = hlog(x)`, then you may call `hexp(Huge, ix)` instead.

If you already know `lx = log(x)`, then you may call `exp(Huge, ix)` instead.
"""
struct Huge{T<:Real} <: Real
    hlog::T
    global hexp(::Type{Huge{T}}, x::T) where {T<:Real} = new{T}(x)
end

const HugeFloat16 = Huge{Float16}
const HugeFloat32 = Huge{Float32}
const HugeFloat64 = Huge{Float64}

decon(x::Huge) = (ix = hlog(x); (signbit(ix), abs(ix)))
recon(sb::Bool, ix::Real) = hexp(Huge, ifelse(sb, -ix, +ix))

"""
    hexp(Huge, x)

Compute `hexp(x)` but the result is stored exactly as a `Huge`.
"""
hexp(::Type{Huge{T}}, x::Real) where {T<:Real} = hexp(Huge{T}, convert(T, x))
hexp(::Type{Huge}, x::T) where {T<:Real} = hexp(Huge{T}, x)

hlog(x::Huge) = x.hlog
hexp(x::Huge) = hexp(Huge, hexp(hlog(x)))

Base.convert(::Type{Huge{T}}, x::Huge{T}) where {T} = x
Base.convert(::Type{Huge{T}}, x::Real) where {T} = hexp(Huge{T}, hlog(x))

Base.convert(::Type{Huge}, x::Huge) = x
Base.convert(::Type{Huge}, x::Real) = hexp(Huge, hlog(x))

Huge(x::Real) = convert(Huge, x)

Huge{T}(x::Real) where {T} = convert(Huge{T}, x)

Base.AbstractFloat(x::Huge) = AbstractFloat(hexp(AbstractFloat(hlog(x))))

Base.BigFloat(x::Huge) = BigFloat(hexp(BigFloat(hlog(x))))

Base.promote_rule(::Type{Huge{T}}, ::Type{Huge{S}}) where {T<:Real, S<:Real} = Huge{promote_type(T, S)}
Base.promote_rule(::Type{Huge{T}}, ::Type{S}) where {T<:Real, S<:Real} = promote_type(Huge{T}, typeof(Huge(zero(S))))


### COMPARISONS / PREDICATES

for op in [:signbit, :isfinite, :isinf, :isnan, :iszero, :isone]
    @eval Base.$op(x::Huge) = $op(hlog(x))
end

for op in [:isequal, :isless, :cmp, :(==), :(!=), :(<), :(<=), :(>), :(>=)]
    @eval Base.$op(x::Huge, y::Huge) = $op(hlog(x), hlog(y))
end

for op in [:sign, :abs, :(+), :(-), :inv, :nextfloat, :prevfloat]
    @eval Base.$op(x::Huge) = hexp(Huge, $op(hlog(x)))
end

for op in [:zero, :one, :typemin, :typemax]
    @eval Base.$op(::Type{Huge{T}}) where {T} = hexp(Huge, $op(T))
end

for op in [:nextfloat, :prevfloat]
    @eval Base.$op(x::Huge, n::Integer) = hexp(Huge, $op(hlog(x), n))
end

Base.zero(::Type{Huge}) = zero(Huge{Int})

Base.one(::Type{Huge}) = one(Huge{Int})

# We hash the inner value, which for free means that hash(Huge(x)) == hash(x) for x in
# [-Inf, -1, 0, 1, Inf], which are also the only values likely to overlap with other types.
Base.hash(x::Huge, h::UInt) = hash(hlog(x), h)


### TYPES

Base.big(::Type{Huge}) = Huge{BigFloat}
Base.big(::Type{Huge{T}}) where {T} = Huge{big(T)}
Base.big(x::Huge) = convert(big(typeof(x)), x)

Base.widen(::Type{Huge{T}}) where {T} = Huge{widen(T)}
Base.widen(x::Huge) = convert(widen(typeof(x)), x)


### ARITHMETIC

function Base.:(*)(x::Huge{T}, y::Huge{T}) where {T}
    sx, ix = decon(x)
    sy, iy = decon(y)
    sr = xor(sx, sy)
    ir = invh(h(ix) + h(iy))
    return recon(sr, ir)
end

function Base.:(/)(x::Huge{T}, y::Huge{T}) where {T}
    sx, ix = decon(x)
    sy, iy = decon(y)
    sr = xor(sx, sy)
    ir = invh(h(ix) - h(iy))
    return recon(sr, ir)
end

function Base.:(^)(x::Huge, p::Real) where {T}
    sx, ix = decon(x)
    if sx && !iszero(ix) && !isnan(ix)
        if isinteger(p)
            sr = isodd(p)
        else
            throw(DomainError(x, "Can only take integer powers of negative Huge numbers."))
        end
    else
        sr = false
    end
    ir = invh(h(ix) * p)
    return recon(sr, ir)
end

function Base.:(^)(x::Huge, p::Integer) where {T}
    sx, ix = decon(x)
    sr = sx && isodd(p)
    ir = invh(h(ix) * p)
    return recon(sr, ir)    
end

Base.:(+)(x::Huge{T}, y::Huge{T}) where {T} = _add(decon(x)..., decon(y)...)

function _add(sx, ix, sy, iy)
    se = xor(sx, sy)
    sr = se ? xor(ix > iy, sy) : sx
    a, b = minmax(ix, iy)
    ha = h(a)
    hb = h(b)
    c = ha - hb
    if se
        if c < oftype(c, -1)
            d = log(-expm1(c))
        else
            d = log1p(-exp(c))
        end
    else
        d = log1p(exp(c))
    end
    ir = invh(hb + d)
    return recon(sr, ir)    
end

Base.:(-)(x::Huge{T}, y::Huge{T}) where {T} = _sub(decon(x)..., decon(y)...)

_sub(sx, ix, sy, iy) = _add(sx, ix, !sy, iy)

function Base.log(x::Huge)
    ix = hlog(x)
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
    return logx / log(oftype(logx, 10))
end

function Base.exp(x::Huge)
    return hexp(Huge, invh(hexp(hlog(x))))
end

function Base.exp(::Type{Huge}, x::Real)
    if signbit(x)
        ir = inv(one(x) - x)
    else
        ir = x + inv(inv(one(x)))
    end
    return hexp(Huge, ir)
end

function Base.exp(::Type{Huge{T}}, x::Real) where {T}
    if signbit(x)
        ir = inv(one(T) - convert(T, x))
    else
        ir = convert(T, x) + one(T)
    end
    return hexp(Huge{T}, ir)
end


## IO

function Base.show(io::IO, x::Huge)
    print(io, "hexp(")
    if get(io, :typeinfo, Any) != typeof(x)
        show(io, typeof(x))
        print(io, ", ")
    end
    show(io, hlog(x))
    print(io, ")")
end

function Base.write(io::IO, x::Huge)
    write(io, hlog(x))
end

function Base.read(io::IO, ::Type{Huge{T}}) where {T}
    hexp(Huge{T}, read(io, T))
end

end # module
