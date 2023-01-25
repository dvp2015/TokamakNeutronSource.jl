### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ 26433bb2-6ff1-11ed-0943-01e79517b4f5
begin
	using Pkg
	Pkg.activate(temp=true)
	packages = [
		"Cuba",
		"GLMakie",
	]
	Pkg.add(packages)
	devdir = Pkg.devdir()
	eqdsk_pkg_dir = joinpath(devdir, "EQDSKReader.jl")
	tns_pkg_dir = joinpath(devdir, "TokamakNeutronSource")
	Pkg.add(path=eqdsk_pkg_dir)
	Pkg.add(path=tns_pkg_dir)
	using EQDSKReader, TokamakNeutronSource
end

# ╔═╡ 9597cccc-5f54-4260-bf2e-2fcf4f3a169a
using Cuba, GLMakie

# ╔═╡ e8552241-d48d-4469-abae-484db725af36
using TokamakNeutronSource.PlasmaDistributions

# ╔═╡ 5f37ce53-477f-498d-b1ea-613fae26f3f8
begin
	eqdsk_path = joinpath(tns_pkg_dir, "test", "data", "beforeTQ.eqdsk")
	excel_path = joinpath(tns_pkg_dir, "test", "data", "TRT_215_8T_NBI.xlsx")
	@assert all(isfile, [eqdsk_path, excel_path])
end

# ╔═╡ 50b9960f-7673-4975-b24b-d96ebb16d50f
eqdsk = Content(eqdsk_path)

# ╔═╡ 7e316d24-32ca-4084-b7ec-fa0b3cee1894
excel = load_excel(excel_path)

# ╔═╡ 11d8369b-a698-4cdb-96ed-5db4244f099e
distr = DDDistribution(eqdsk, excel)

# ╔═╡ ad233063-b49b-49f2-a6d9-511324b30486
let
	using TokamakNeutronSource.Integrations
	total = total_yield(distr)[1]
	rmin, rmax, zmin, zmax = domain(distr)
	moment0 = torroidal_segment_moment_0((r,z) -> I(distr, r, z), (rmin, rmax), (zmin, zmax))[1]
	@assert moment0 ≈ total
end

# ╔═╡ 0f292ede-4a1f-4387-9f49-e0ac8ee993de
begin
	r = rpoints(eqdsk)
	z = zpoints(eqdsk)
	_I = I(distr, r, z)
	extrema(_I)
end

# ╔═╡ 5d1b8d4f-bccd-43c1-81ff-4c186b3b5367
typeof(I)

# ╔═╡ e7b4f8b3-f3f0-4d5c-9080-eea93e01dd17
function plot_neutron_source(eqdsk::Content, distr::AbstractDistribution)
	f = Figure(resolution=(600, 800))
	ax = Axis(
		f[1,1]; 
		xlabel=L"R,m", 
		ylabel=L"Z,m", 
		aspect=DataAspect(),
		title=L"$$DD neutron source intensity"
	)
	r = rpoints(eqdsk)
	z = zpoints(eqdsk)
	xlims!(ax, r[1], r[end])
	ylims!(ax, z[1], z[end])
	_I = I(distr, r, z)
	levels=range(0, 1.3e10, length=14)
	cntr = contourf!(ax,
		r, z, _I,
		levels=levels, 
		colormap=:bamako, 
		linestyle="-",
		linecolor=:black,
		linewidth=2
	)
	rbbs_points = [Point2f(x,y) for (x,y) in zip(eqdsk.rbbbs, eqdsk.zbbbs)]
	rlim_points = [Point2f(x,y) for (x,y) in zip(eqdsk.rlim, eqdsk.zlim)]
	scatter!(ax, eqdsk.rmaxis, eqdsk.zmaxis, color=:gray60, marker=:cross, label="Magnetic axis")
	lines!(ax, rbbs_points, label="Plasma boundary")
	lines!(ax, rlim_points, label="Limiter")
	cb=Colorbar(f[1,2], cntr, label=L"$I_{DD}(R,Z)$, $cm^{-3}s^{-1}", ticks=levels[1:2:end])
	cb.tellheight = true
	cb.tellwidth = true
	axislegend(ax)
	# axs = Axis3(f[1,2])
	# surface!(axs, r,z,I)
	f	
end

# ╔═╡ 6888f022-b9dc-43f6-8759-772f97afc734
plot_neutron_source(eqdsk, distr)

# ╔═╡ ae39aefc-80ee-437d-9bb6-e763206e4339
md"""
# Define appropiate mesh

The mesh quality criteria - the variance of I(r,z) in a cell not more than allowed_variance.
"""


# ╔═╡ cbb99594-30c3-4ccf-9251-677fc6755750
function cell_relative_variance(A, i, j)
	subA = @view A[i:i+1,j:j+1]
	a, b = extrema(subA)
	nom = b - a	
	denom = sum(subA)
	denom <= 0.0 && return 1.0
	res = 4nom / denom
	min(1.0, res)
end

# ╔═╡ f453daec-8564-48e4-8869-834d3ef64fd7
function variance_on_mesh(f, rbins, zbins)
	nr, nz = length(rbins), length(zbins)
	matrix = zeros(nr, nz)
	for i in 1:nr-1, j in 1:nz-1
		matrix[i,j] = integrate_torroidal_segment(f, rbins[i:i+1], zbins[j:j+1])[1] 
	end
	[cell_relative_variance(matrix, i, j) for i = 1:size(matrix, 1)-1, j = 1:size(matrix,2)-1 ]
end

# ╔═╡ 643e388a-18b3-4c16-bc10-d6f0e455fafd
let
	_I(r,z) = I(distr, r, z)
	rmin, rmax, zmin, zmax = domain(distr)
	factor(n::Int) = n  # Int(floor(0.7n))
	r  = range(rmin, rmax, length=factor(65))
	z  = range(zmin, zmax, length=factor(129))
	# map(length, [r, z])
	rmids = 0.5(r[1:end-1] .+ r[2:end])
	zmids = 0.5(z[1:end-1] .+ z[2:end])
	vom = variance_on_mesh(_I, r, z)
	f = Figure(resolution=(600, 800))
	ax = Axis(
		f[1,1]; 
		xlabel=L"R,m", 
		ylabel=L"Z,m", 
		aspect=DataAspect(),
		title=L"$$DD neutron source intensity variance"
	)
	xlims!(ax, rmin, rmax)
	ylims!(ax, zmin, zmax)
	# top = ceil(100 * maximum(vom))
	# levels=range(0, top*0.01, length=Int(top)+1)
	# levels=range(0, 1.3, length=11)
	# levels = range(0,0.3, length=11)
	cntr = contourf!(ax,
		rmids, zmids, vom,
		# levels=levels, 
		colormap=:bamako, 
		linestyle="-",
		linecolor=:black,
		linewidth=2
	)
	rbbs_points = [Point2f(x,y) for (x,y) in zip(eqdsk.rbbbs, eqdsk.zbbbs)]
	rlim_points = [Point2f(x,y) for (x,y) in zip(eqdsk.rlim, eqdsk.zlim)]
	scatter!(ax, eqdsk.rmaxis, eqdsk.zmaxis, color=:gray60, marker=:cross, label="Magnetic axis")
	lines!(ax, rbbs_points, label="Plasma boundary")
	lines!(ax, rlim_points, label="Limiter")
	cb=Colorbar(
		f[1,2], 
		cntr, 
		label=L"Variance", 
		# ticks=levels[2:2:end]
	)
	cb.tellheight = true
	cb.tellwidth = true
	axislegend(ax)
	# axs = Axis3(f[1,2])
	# surface!(axs, r,z,I)
	f	
end

# ╔═╡ e453e4a7-7fdb-4adb-9b31-d1f89f83bf13
md"""
Total yield and moment-0 of a distribution are the same.
"""

# ╔═╡ f328d540-9e59-47dd-8169-9b022be9e7f6
let
	rmin, rmax, zmin, zmax = domain(distr)
	torroidal_segment_moment_0((r,z) -> I(distr, r, z), (rmin, rmax), (zmin, zmax))[1]
end

# ╔═╡ e2dde790-8e0f-4189-9971-a366ff2d6526
md"""
	Moment 1 is close magnetic axes coordinates.
"""

# ╔═╡ e8381d26-3899-4f63-a905-aa7a47333401
let
	rmin, rmax, zmin, zmax = domain(distr)
	torroidal_segment_moment_1((r,z) -> I(distr, r, z), (rmin, rmax), (zmin, zmax))[1]
end

# ╔═╡ a952abea-f29e-46bf-9e7e-8e9bd16b1b2c
eqdsk.rmaxis, eqdsk.zmaxis

# ╔═╡ 6eec651c-e7ee-40e6-89f7-0b37af8e6c15
md"""
## Compute params for SDEF computation.
"""

# ╔═╡ 187df34e-77a7-4c0f-9d10-ba1edfbdd41f
let
	rmin, rmax, zmin, zmax = domain(distr)
	r  = range(rmin, rmax, length=65)
	z  = range(zmin, zmax, length=129)
	100 .* extrema(diff(r)), 100 .* extrema(diff(z))
end

# ╔═╡ 1f5d7e83-7a13-4046-a7cc-7aa03e293d18
function compute_sdef_values(distr, nr, nz)
	rmin, rmax, zmin, zmax = domain(distr)
	f(r,z) = I(distr, r, z)
	rbins  = range(rmin, rmax, length=nr)
	zbins  = range(zmin, zmax, length=nz)
	src = zeros(nr, nz)
	for i in 1:nr-1, j in 1:nz-1
		src[i,j] = integrate_torroidal_segment(f, rbins[i:i+1], zbins[j:j+1])[1] 
	end
	rbins, zbins, src ./ sum(src)
end

# ╔═╡ 216554fb-34d9-4ab8-83e3-bab6f438c5df
rbins, zbins, src = compute_sdef_values(distr, eqdsk.nw, eqdsk.nh);

# ╔═╡ 6d363253-a395-4f6b-91e9-21e8f9fcabb8
maximum(src), sum(src), size(src), length(src)

# ╔═╡ b2df9bd8-64d6-4e86-9cbd-2b53e09cdc87
let
	ci = argmax(src)
	rbins[ci[1]], zbins[ci[2]]
end

# ╔═╡ Cell order:
# ╠═26433bb2-6ff1-11ed-0943-01e79517b4f5
# ╠═9597cccc-5f54-4260-bf2e-2fcf4f3a169a
# ╠═5f37ce53-477f-498d-b1ea-613fae26f3f8
# ╠═50b9960f-7673-4975-b24b-d96ebb16d50f
# ╠═e8552241-d48d-4469-abae-484db725af36
# ╠═7e316d24-32ca-4084-b7ec-fa0b3cee1894
# ╠═11d8369b-a698-4cdb-96ed-5db4244f099e
# ╠═0f292ede-4a1f-4387-9f49-e0ac8ee993de
# ╠═5d1b8d4f-bccd-43c1-81ff-4c186b3b5367
# ╠═e7b4f8b3-f3f0-4d5c-9080-eea93e01dd17
# ╠═6888f022-b9dc-43f6-8759-772f97afc734
# ╠═ae39aefc-80ee-437d-9bb6-e763206e4339
# ╠═f453daec-8564-48e4-8869-834d3ef64fd7
# ╠═643e388a-18b3-4c16-bc10-d6f0e455fafd
# ╟─cbb99594-30c3-4ccf-9251-677fc6755750
# ╟─e453e4a7-7fdb-4adb-9b31-d1f89f83bf13
# ╟─ad233063-b49b-49f2-a6d9-511324b30486
# ╠═f328d540-9e59-47dd-8169-9b022be9e7f6
# ╟─e2dde790-8e0f-4189-9971-a366ff2d6526
# ╠═e8381d26-3899-4f63-a905-aa7a47333401
# ╠═a952abea-f29e-46bf-9e7e-8e9bd16b1b2c
# ╠═6eec651c-e7ee-40e6-89f7-0b37af8e6c15
# ╠═187df34e-77a7-4c0f-9d10-ba1edfbdd41f
# ╠═1f5d7e83-7a13-4046-a7cc-7aa03e293d18
# ╠═216554fb-34d9-4ab8-83e3-bab6f438c5df
# ╠═6d363253-a395-4f6b-91e9-21e8f9fcabb8
# ╠═b2df9bd8-64d6-4e86-9cbd-2b53e09cdc87
