#   II. Cuspidality analysis: the cuspidal case
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡
# 
#   In this session, we will analyse the cuspidality of elementary robots using
#   AlgebraicSolving.jl.

using Pkg, Plots, AlgebraicSolving
import AlgebraicSolving.change_ringvar, AlgebraicSolving.gens, AlgebraicSolving.factor

#   A 3 revolute joint orthogonal robot:
#   ––––––––––––––––––––––––––––––––––––
# 
#   (Image: A 3R robot)

#   We start with a robot whose algebraic kinematic map is: :$
#   \begin{array}{cccc} \mathcal{K} \colon& V &\longrightarrow & \mathbb{R}^3 \
#   &(c1,c2, c3, s1, s2, s3) &\longmapsto &\big(z1(\bm{c},\bm{s}),
#   z2(\bm{c},\bm{s}),z_3(\bm{c},\bm{s})\big) \end{array} :$

#   where \left\{\begin{array}{l}  z_1 = \frac{1}{2}c_1c_2(3c_3 + 4) -
#   \frac{1}{2}s_1(3s_3 + 2) + c_1\\[.5em]  z_2 = \frac{1}{2}s_1c_2(3c_3 + 4) +
#   \frac{1}{2}c_1(3s_3 + 2) + s_1\\[[.5em] z_3 = -\frac{1}{2}s_2(3c_3 + 4) 
#   \end{array}\right.

#   A. System Definition
#   ––––––––––––––––––––
# 
#   First, initialize the multivariate polynomial ring over the rationals (QQ)
#   for our 6 variables: c1, c2, c3, s1, s2, s3.

R, (c1, c2, c3, s1, s2, s3) = polynomial_ring(QQ, [:c1,:c2,:c3,:s1,:s2,:s3])

#   Now, we define our kinematic map K and the ideal V constraining the
#   (co)sines (c_i^2 + s_i^2 - 1 = 0).

K = 1//2 * [ c1*c2*(3*c3+4) - s1*(3*s3+2) + c1, s1*c2*(3*c3+4) + c1*(3*s3+2) + s1, -s2*(3*c3+4) ]
V = Ideal([c1^2+s1^2-1, c2^2+s2^2-1, c3^2+s3^2-1])

[K, V]

#   <details> <summary><b>► Show Solution</b></summary>
# 
#   V = Ideal([ c1^2 + s1^2 - 1, c2^2 + s2^2 - 1, c3^2 + s3^2 - 1 ])
# 
#   </details>

dimension(V), hilbert_degree(V)

#   B. Polar Varieties and Critical Points
#   ––––––––––––––––––––––––––––––––––––––
# 
#   Use the computepolar function to extract the critical points with respect to
#   the map K.

# Compute the polar variety
critK = computepolar(1:3, V, phi=K)

#   <details> <summary><b>► Show Solution</b></summary>
# 
#   critK = computepolar(1:3, V, phi=K)
# 
#   </details>

dimension(Ideal(critK)), hilbert_degree(Ideal(critK))

#   Factoring defining polynomials can often lead to (substantial)
#   simplifications

fact = AlgebraicSolving.factor(critK[end])
L_fact = getindex.(collect(fact), 1)

#   Pick the right factor(s)

critK[end] = L_fact[end]

#   <details> <summary><b>► Show Solution</b></summary>
# 
#   critK[end] = L_fact[end]
# 
#   </details>

#   C. The Direct Approach
#   ––––––––––––––––––––––
# 
#   Let's introduce 3 new variables (z_1, z_2, z_3) and attempt standard
#   elimination.
# 
#   We must first inject polynomials in the extended ring.

K_z = [ change_ringvar(f, [:c1,:c2,:c3,:s1,:s2,:s3,:z1,:z2,:z3]) for f in K ]
critK_z = [ change_ringvar(f, [:c1,:c2,:c3,:s1,:s2,:s3,:z1,:z2,:z3]) for f in critK ]

#   We get the new variables

Z = AlgebraicSolving.gens(parent(K_z[1]))[end-2:end]

#   Set Z[1], Z[2], Z[3] to be the image of K on V by adding new equations
#   encoding it

new_eqs = [Z[1] - K_z[1], Z[2] - K_z[2], Z[3] - K_z[3]]
KcritK = Ideal(vcat(critK_z, new_eqs))

#   <details> <summary><b>► Show Solution</b></summary>
# 
#   [Z[1] - K_z[1], Z[2] - K_z[2], Z[3] - K_z[3]]
# 
#   </details>

#   We now project on the z_i by elimination and inject in the polynomial ring
#   only in the z_i

# Eliminate the variable
G = eliminate(KcritK, 6)
svalK = change_ringvar(G)

#   We now attempt to compute one point in each connected component of
#   \mathbb{R}^3 \setminus \mathrm{sval}(K).

Q_sval = points_per_components(
QQMPolyRingElem[],
QQMPolyRingElem[],
svalK,
info_level=1)

#   <details> <summary><b>► Show Solution</b></summary>
# 
#   Q_sval = points_per_components(
#   QQMPolyRingElem[],
#   QQMPolyRingElem[],
#   svalK,
#   info_level=1)
# 
#   </details>

#   Axial symmetry w.r.t. \theta_1

# Configuration space
groebner_basis(Ideal(critK))

# Workspace (cylindrical coordinates)
A, (rho, cp, sp, z) = polynomial_ring(QQ, [:rho,:cp,:sp,:z])
groebner_basis(Ideal([cp^2+sp^2-1, svalK[1](rho*cp, rho*sp, z)]))

#   D. A geometric simplification
#   –––––––––––––––––––––––––––––
# 
#   By intelligently changing our coordinate system (using \rho and z), we can
#   bypass the previous explosion in complexity.

R_polar, (c2,c3,s2,s3,rho,z) = polynomial_ring(QQ, [:c2,:c3,:s2,:s3,:rho,:z])
K_polar, critK_polar = [ [ f(1,c2,c3,0,s2,s3) for f in E if !iszero(f(1,c2,c3,0,s2,s3))] for E in [K, critK] ]

# We now project on (rho, z)
KcritK_polar = Ideal(vcat([z - K_polar[end], rho^2 - K_polar[1]^2 - K_polar[2]^2], critK_polar))

# Eliminate the variables
G = eliminate(KcritK_polar, 4)
svalK_polar = change_ringvar(G)

Wsval_polar = points_per_components(QQMPolyRingElem[], [rho], svalK_polar, info_level=2)

#   E. Roadmap computation
#   ––––––––––––––––––––––
# 
#   We now need to compute a roadmap for the semi-algebraic set defined by V
#   \setminus \mathrm{crit}(\mathcal{K}, V) i.e. by :$ ci^2 + si^2 =1
#   \text{\quad and \quad} \texttt{critK[end]} \neq 0 :$
# 
#   To reduce to algebraic sets, we use here for simplicity the Rabinowitsch's
#   trick to convert disequalities into equalities. We represent the
#   semi-algebraic set as a projection by introducing a new variable u: :$ { f
#   \neq 0 } \longrightarrow {u\cdot f = 1}. :$
# 
#   NB: In practice, this leads to complexity explosion as complexities are
#   intrinsically exponential in a polynomial of number of variables.

R_rab,(u,c2,c3,s2,s3) = polynomial_ring(QQ, [:u,:c2,:c3,:s2,:s3])
critK_rab = [ f(c2,c3,s2,s3,0,0) for f in critK_polar ]
K1_rab = [ f(c2,c3,s2,s3,0,0) for f in K_polar ]

# Lifting of the sa set as an algebraic set
Vmcrit_rab = vcat(critK_rab[1:end-1], 1 - u*critK_rab[end])

#   We recompute the candidate points

# Fiber equations
Wsval_rab = [ vcat(Vmcrit_rab, sum(rfib)//2 - K1_rab[1]^2 - K1_rab[2]^2, sum(zfib)//2 - K1_rab[end]) for (rfib, zfib) in Wsval_polar ]
# Real root isolation
#Wcusp_isolated_rab = filter(!isempty, [ real_solutions(Ideal(w), interval=true)  for w in Wcusp_rab ])

RM_rab = roadmap(Ideal(Vmcrit_rab), Ideal(Wsval_rab[end]), info_level=1)

hilbert_degree.(all_eqs(RM_rab))

GRM_rab = curve_arrangement_graph(all_eqs(RM_rab), Ideal(Wsval_rab[end]), generic=15, v=3)

group_by_component(GRM_rab)

number_of_connected_components(GRM_rab)

#   Bonus: vizualization

function plot_graph(G)
  P = build_graph_data(G)
  plot(legend=false)
  E = P.edge_group
  plot!(E.edges, color=E.color, width = E.width)
  [scatter!(V.vertices, color=V.color, marker = V.marker) for V in P.point_groups]
  plot!()
  #gui()
end

function plot_graphs(CG)
  CP = build_graph_data(CG, vert=false)
  plot(legend=false)
  for P in CP
      E = P.edge_group
      plot!(E.edges, color=E.color, width = E.width)
      [ scatter!(V.vertices, color=V.color, marker = V.marker) for V in P.point_groups ]
  end
  plot!()
  #gui()
end
plot_graphs(connected_components(GRM_rab))