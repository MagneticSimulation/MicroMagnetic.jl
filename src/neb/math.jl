function spherical2cartesian(spherical, n_total)
    dev = get_backend(spherical)
    cartesian = KernelAbstractions.zeros(dev, Float[], 3*n_total)
    spherical2cartesian_kernel!(dev, groupsize[])(spherical, cartesian, ndrange=n_total)
    KernelAbstractions.synchronize(dev)
    return cartesian
end

function cartesian2spherical(cartesian, n_total)
    dev = get_backend(cartesian)
    spherical = KernelAbstractions.zeros(dev, Float[], 3*n_total)
    cartesian2spherical_kernel!(dev, groupsize[])(cartesian, spherical, ndrange=n_total)
    KernelAbstractions.synchronize(dev)
    return spherical
end

function inner_product(x1, x2, n_total)
    dev = get_backend(x1)
    res = KernelAbstractions.zeros(dev, Float[], n_total)
    inner_product_kernel!(dev, groupsize[])(x1, x2, res, ndrange=n_total)
    KernelAbstractions.synchronize(dev)
    return res
end

function slerp(x1, x2, t, n_total::Int)
    dev = get_backend(x1)
    res = KernelAbstractions.zeros(dev, Float[], 3*n_total)
    slerp_kernel!(dev, groupsize[])(x1, x2, t, res, ndrange=n_total)
    KernelAbstractions.synchronize(dev)
    return res
end
