export GetGzBethe
export GetGzBetheDim
export GetGzBetheUniformScaling
export GetGzBetheFromSelf


@inline GetGiwBethe( w ; D=1 ) = 2.0/D/D * ( w - sign(imag(w)) * sqrt(Complex( w*w - D*D )) )
@inline GetGwBethe(  w ; D=1 ) = 2.0/D/D * ( w - sign(real(w)) * sqrt(Complex( w*w - D*D )) )

@inline GetGzBethe(  w ; D=1 ) = 2.0*(w + sqrt(Complex(1-w^2))*(log(Complex(1-w)) - log(Complex(-1+w)))/pi)

@inline function GetGzBetheDim(  w , dim ; D=1 )
    return collect(I(dim)) * 2.0*( w + sqrt(Complex(1-w^2))*(log(Complex(1-w)) - log(Complex(-1+w)))/pi)
end

@inline function GetGzBetheUniformScaling(  w ; D=1 )
    return I*2.0*( w + sqrt(Complex(1-w^2))*(log(Complex(1-w)) - log(Complex(-1+w)))/pi)
end

@inline function GetGzBetheUniformScalingDim(  w, dim ; D=1 )
    return I(dim)*2.0*( w + sqrt(Complex(1-w^2))*(log(Complex(1-w)) - log(Complex(-1+w)))/pi)
end

@inline function GetGzBetheFromSelf(  w::ComplexF64 , self ; D=1 )
    z   = w - self
    return GetGzBethe( z ; D=D)
end

@inline function GetGzBetheFromSelf(  w::Vector{ComplexF64} , self_w ; D=1 )
    nfreq   = length(w)
    return [ GetGzBetheFromSelf( w[ifreq] , self_w[ifreq] ; D=D) for ifreq in 1:nfreq ]
end
