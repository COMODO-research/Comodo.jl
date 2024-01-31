using Rotations

mutable struct braccio
    α₁ :: Float64
    α₂ :: Float64
    α₃ :: Float64
    α₄ :: Float64
    α₅ :: Float64
    α₆ :: Float64
    P :: Vector{Float64}
    Q :: RotXYZ{Float64}
end

braccio(α₁,α₂,α₃,α₄,α₅,α₆,P,Q)=braccio(α₁=0.0,α₂=0.0,α₃=0.0,α₄=0.0,α₅=0.0,α₆=0.0,P=[0.0,0.0,0.0],Q = RotXYZ(0.0,0.0,0.0)) 