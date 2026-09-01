using Pkg, Plots, AlgebraicSolving
import AlgebraicSolving.change_ringvar, AlgebraicSolving.gens, AlgebraicSolving.factor

#   I. Warm-up
#   ≡≡≡≡≡≡≡≡≡≡

#   1.Algebraic geometry
#   ====================

# Define the ring of polynomials in (x1 , x2, x3), with rational coefficients
R, (x1, x2, x3) = polynomial_ring(QQ,[:x1,:x2,:x3])

# Ideal of definition of the 3-sphere
I = Ideal([x1^2 + x2^2 + x3^2 - 1])

dimension(I), hilbert_degree(I)

# Generic fiber of the sphere
J = Ideal(vcat(I.gens, [x1, x3]))

# Compute the finite (?) points
P = rational_parametrization(J)

# Project on (x2,x3)
projI = eliminate(I, 1)

#   2. Real algebraic geometry: projection
#   ======================================

#   Compute the <b>real</b> projection of the sphere on (x_2,x_3)
#   -------------------------------------------------------------

# Set x2,x3 as a parameters: define QQ(x2,x3)[x1]
A, (x2_p, x3_p) = polynomial_ring(QQ, [:x2,:x3])
K = fraction_field(A)
B, (x1_p,) = polynomial_ring(K, [:x1])

# define parametric Ideal
I_p = ParametricIdeal([I.gens[1](x1_p,x2_p,x3_p)])

RC_I = real_root_classification(I_p, [B(1)], info_level=5)

RC_I.S

#   Projection of the sphere on x_3: polar varieties
#   ------------------------------------------------

# Set x3 as a parameter: define QQ(x3)[x1,x2]
A, (x3_p, ) = polynomial_ring(QQ, [:x3])
K = fraction_field(A)
B, (x1_p,x2_p) = polynomial_ring(K, [:x1,:x2])
I_p = ParametricIdeal([I.gens[1](x1_p,x2_p,x3_p)])

RC_I = real_root_classification(I_p, [B(1)])

#   Polar variety of the projection on x_1.
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅
# 
#   In a fiber \longrightarrow is finite and reaches every connected component.

# For type reasons we must pass by the generators and provide the dimension in Q(x3)[x1,x2]
Itemp = Ideal(I_p.gens) # classical ideal
Itemp.dim = 1
W = computepolar(1:1, Itemp)

I_W = ParametricIdeal(W)
RC_IW = real_root_classification(I_W, [B(1)], info_level=5)

RC_IW.S

#   Non-empty for 1-x_3^2 > 0 i.e. -1 < x_3 < 1

#   3. Sampling
#   ===========

#   Proof of Turan's inequality for Legendre polynomials
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅
# 
#   $ Pn(x)^2 - P{n-1}(x)P_{n+1}(x) \geq 0$ for x \in [-1,1], \; n\geq 1
# 
#   with (n+2)P_{n+2} - (2n+3)xP_{n+1} + (n+1)P_{n} = 0

R,(n,x,p0,p1) = polynomial_ring(QQ, [:n,:x,:p0,:p1]);

# Parameters non-negative constraints
params = [n-1, x+1, 1-x];

# Recurrence for substitution
npm1, n2p2 = -(n+1)*p1+(2*n+1)*x*p0, (2*n+3)*x*p1-(n+1)*p0;

# Turan's non-negative conditions
Tu0 = n*p0^2-npm1*p1;
Tu1 = (n+2)*p1^2 - p0*n2p2;

X = points_per_components(QQMPolyRingElem[], vcat(params, [Tu0, -Tu1]), [Tu1], info_level=5)

#   Invariant generation for polynomial loops
#   ⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅

#   (x_1,x_2,x_3)=(1,-1,1)
# 
#   \textbf{while } true \textbf{ do}
# 
#   \qquad\begin{pmatrix}
# x_1 \\ x_2 \\ x_3
# \end{pmatrix}
# \leftarrow
# \begin{pmatrix}
# 2x_1 \\ -3x_2 \\ 4x_3
# \end{pmatrix}

R,(x1,x2,x3) = polynomial_ring(QQ, [:x1,:x2,:x3])

#   Non inductivity of the invariant 2x_2+4x_3-1>0

Init = map(QQ,[1,-1,1])
Inv = 2*x2+4*x3-1
F = [2*x1, -3*x2, 4*x3]
InvF = Inv(F...)

println("Initialization: ", Inv(Init...) > 0)
println("Induction: ", isempty(points_per_components(QQMPolyRingElem[], [Inv, -InvF], [Inv, InvF])))

#   But 2-inductivity

println("Initializations: ", Inv(Init...) > 0, ", ", InvF(Init...) > 0)
InvFF = InvF(F...)
println("Induction: ", isempty(points_per_components(QQMPolyRingElem[], [Inv, InvF, -InvFF], [Inv, InvF, InvFF])))

A,(c0,c1,c2) = polynomial_ring(QQ,[:c0,:c1,:c2])
K = fraction_field(A)
R, (x1,x2) = polynomial_ring(K, [:x1,:x2])

Init = QQ.([1,-1])
Inv = c0 + c1*x1 + c2*x2
F = [10*x1-8*x2, 6*x1-4*x2]
InvF = Inv(F...)

I_inv = ParametricIdeal([Inv])

lol = Ideal(I_inv.gens)
lol.dim = 1
J_inv = ParametricIdeal(computepolar(1:1, lol, phi=[x1^2 + x2^2]))

RR_J = real_root_classification(J_inv, [InvF], info_level=5)

S1, S2 = RR_J.S

#   The inequations

S1.ineqs

#   The polynomial whose sign must be fixed so that there is no root above

S2.polys[1], first.(S2.signs[map(==(0), S2.counts)])

S1_spec, S2_spec = [ [ f(c2-c1,c1,c2) for f in E] for E in [S1.ineqs, S2.polys] ]

P_inv = points_per_components(QQMPolyRingElem[], S2_spec, S1_spec, info_level=5)

unique([ sum.([c2-c1, c1, c2]) for (c1,c2) in P_inv ])

#   4.Algebraic curves
#   ==================

function plot_graph(P)
         plot(legend=false)
         E = P.edge_group
         plot!(E.edges, color=E.color, width = E.width)
         [scatter!(V.vertices, color=V.color, marker = V.marker) for V in P.point_groups]
         plot!()
end

function plot_graphs(CP)
         plot(legend=false)
         for P in CP
             E = P.edge_group
             plot!(E.edges, color=E.color, width = E.width)
             [ scatter!(V.vertices, color=V.color, marker = V.marker) for V in P.point_groups ]
         end
         plot!()
end

R, (x1,x2,x3) = polynomial_ring(QQ, [:x1, :x2,:x3])
IC = Ideal([x2 - x1^2 + 1, x3 - x1*x2])

G = curve_graph( IC, generic=10, v=5)

plot_graphs(build_graph_data(connected_components(G)))

R, (t,x1,x2,x3) = polynomial_ring(QQ, [:t,:x1,:x2,:x3])
C1 = [x2^2 + x3^2 - 1, x1]
C2 = [x2-x3, x1-1]
C_union = AlgebraicSolving.change_ringvar(eliminate(Ideal(vcat(t.*C1, (1-t).*C2)), 1))

G = curve_graph(Ideal(C_union), generic=20, v=5)

plot_graphs(build_graph_data(connected_components(G)))

G1 = curve_arrangement_graph([Ideal(change_ringvar(C1)), Ideal(change_ringvar(C2))], generic=10)

plot_graphs(build_graph_data(connected_components(G1)))

#   5.Roadmaps
#   ==========

R, (x1,x2,x3) = polynomial_ring(QQ, [:x1,:x2,:x3])

r1, r2 = 3, 1
h = (x1^2+x2^2+x3^2+r1^2-r2^2)^2-4*r1^2*(x1^2+x2^2) + 1//10
I = Ideal([h])

RM = roadmap(I, Ideal([R(1)]), info_level=3)

RM_eq = all_eqs(RM)

LG = curve_arrangement_graph(RM_eq, v=4, generic=5)

plot_graphs(build_graph_data(connected_components(LG)))