#   III. Cuspidality analysis: the non cuspidal case
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡
# 
#   In this session, we will analyse the cuspidality of elementary robots using
#   AlgebraicSolving.jl.

using Pkg, Plots, AlgebraicSolving
import AlgebraicSolving.change_ringvar, AlgebraicSolving.gens, AlgebraicSolving.factor

function log(msg)
    println(msg)
    flush(stdout)
end

#   A 3 revolute joint orthogonal robot:
#   ––––––––––––––––––––––––––––––––––––
# 
#   (Image: A 3R robot)

#   We now consider a robot whose algebraic kinematic map is: :$
#   \begin{array}{cccc} \mathcal{K} \colon& V &\longrightarrow & \mathbb{R}^3 \
#   &(c1,c2, c3, s1, s2, s3) &\longmapsto &\big(z1(\bm{c},\bm{s}),
#   z2(\bm{c},\bm{s}),z_3(\bm{c},\bm{s})\big) \end{array} :$

#   where \left\{\begin{array}{l}  z_1 = \frac{1}{10}c_1c_2(15c_3 + 11) -
#   \frac{1}{10}s_1(15s_3 + 13) + 3c_1\\[.5em]  z_2 = \frac{1}{10}c_1c_2(15c_3 +
#   11) + \frac{1}{10}c_1(15s_3 + 13) + 3s_1\\[[.5em] z_3 =
#   -\frac{1}{10}s_2(15c_3 + 11)  \end{array}\right.

# Polar variety
R,(c1,c2,c3,s1,s2,s3) = polynomial_ring(QQ, [:c1,:c2,:c3,:s1,:s2,:s3])
V = Ideal([ c1^2 + s1^2 - 1, c2^2 + s2^2 - 1, c3^2 + s3^2 - 1 ])
K = 1//10 * [ c1*c2*(15*c3+11) - s1*(15*s3+13) + 3c1, c1*c2*(15*c3+11) - c1*(15*s3+13) + 3s1, -s2*(15*c3+11) ]

#   Compute one point in each connected component of the complement of
#   \mathrm{sval}(V)
#   ----------------------------

#   <details> <summary><b>► Show Solution</b></summary>
# 
#   critK = computepolar(1:3, V, phi=K)
#   
#   R_polar, (c2,c3,s2,s3,rho,z) = polynomial_ring(QQ, [:c2,:c3,:s2,:s3,:rho,:z])
#   K_polar, critK_polar = [ [ f(1,c2,c3,0,s2,s3) for f in E if !iszero(f(1,c2,c3,0,s2,s3))] for E in [K, critK] ]
#   KcritK_polar = Ideal(vcat([z - K_polar[end], rho^2 - K_polar[1]^2 - K_polar[2]^2], critK_polar))
#   
#   G = eliminate(KcritK_polar, 4)
#   svalK_polar = change_ringvar(G)
#   
#   Wsval_polar = points_per_components(QQMPolyRingElem[], [rho], svalK_polar, info_level=2)
# 
#   </details>

#   Compute a roadmap for V \setminus \mathrm{crit}(\mathcal{K}, V)
#   containing the preimage of these above points.
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

#   <details> <summary><b>► Show Solution</b></summary>
# 
#   R_rab,(u,c2,c3,s2,s3) = polynomial_ring(QQ, [:u,:c2,:c3,:s2,:s3])
#   critK_rab = [ f(c2,c3,s2,s3,0,0) for f in critK_polar ]
#   K1_rab = [ f(c2,c3,s2,s3,0,0) for f in K_polar ]
#   
#   Vmcrit_rab = vcat(critK_rab[1:end-1], 1 - u*critK_rab[end])
#   
#   Wcusp_rab = [ vcat(Vmcrit_rab, sum(rfib)//2 - K1_rab[1]^2 - K1_rab[2]^2, sum(zfib)//2 - K1_rab[end]) for (rfib, zfib) in Wsval_polar ]
#   
#   RM_rab = roadmap(Ideal(Vmcrit_rab), Wcusp_isolated_rab, 1)
# 
#   </details>

hilbert_degree.(all_eqs(RM_rab))

#   Compute a homeomorphic graph of the roadmap and conclude on the
#   cuspidality
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅
# 
#   <b>Caution:</b> takes ~1h20, ~1h just for the first curve

#   <details> <summary><b>► Show Solution</b></summary>
# 
#   GRM_rab = curve_arrangement_graph(all_eqs(RM_rab), Ideal.(Wcusp_rab), generic=10, v=2)
# 
#   </details>