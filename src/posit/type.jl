struct Posit{ps, es} <: Real
    s:: Int
    sn:: Int
    k:: Int
    rs:: Int
    ers:: Int
    e:: Int
    fi:: String
    fs:: Int
    f:: Real
    useed:: Real
end

ps(::Posit{p,e}) where {p,e} = p
es(::Posit{p,e}) where {p,e} = e

# function Base.iterate(c::Posit, state = 0)
#     state >= nfields(c) && return
#     return Base.getfield(c, state+1), state+1
# end