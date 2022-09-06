

"""
Return the dimension of an `Interpolation`
"""
@inline getdim(ip::Interpolation{dim}) where {dim} = dim

"""
Return the reference shape of an `Interpolation`
"""
@inline getrefshape(ip::Interpolation{dim,shape}) where {dim,shape} = shape

"""
Return the polynomial order of the `Interpolation`
"""
@inline getorder(ip::Interpolation{dim,shape,order}) where {dim,shape,order} = order

# ------------------------------------
# 拉格朗日多项式插值
include("lagrange.jl")

# ------------------------------------
# 实时存储单元插值系数
include("cellvalues.jl")