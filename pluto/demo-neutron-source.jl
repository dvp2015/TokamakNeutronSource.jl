### A Pluto.jl notebook ###
# v0.19.18

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
	total_yield(distr)[1]
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


# ╔═╡ a1070c6d-3bb2-4cd8-9c49-fc2a079814c2
function cell_variance(a, i, j)
	a, b = extrema(@view a[i:i+1,j:j+1])
	b - a
end

# ╔═╡ cbb99594-30c3-4ccf-9251-677fc6755750
function cell_relative_variance(A, i, j)
	subA = @view A[i:i+1,j:j+1]
	a, b = extrema(subA)
	nom = b - a	
	denom = sum(subA)
	denom <= 0.0 && return 1.0
	res = nom / denom
	1.0 <= res && return 1.0
	4*res
end

# ╔═╡ f453daec-8564-48e4-8869-834d3ef64fd7
function variance_on_mesh(f, rbins, zbins)
	nr, nz = length(rbins), length(zbins)
	matrix = zeros(nr, nz)
	total = torroidal_segment_yield(
		f, [rbins[1],rbins[end]], [zbins[1],zbins[end]]
	)[1]
	for i in 1:nr, j in 1:nz
		matrix[i,j] = f(rbins[i], zbins[j]) 
		# matrix[i,j] = torroidal_segment_yield(f, rbins[i:i+1], zbins[j:j+1])[1] 
	end
	[cell_relative_variance(matrix, i, j) for i = 1:size(matrix, 1)-1, j = 1:size(matrix,2)-1 ]
end

# ╔═╡ 643e388a-18b3-4c16-bc10-d6f0e455fafd
let
	_I(r,z) = I(distr, r, z)
	rmin, rmax, zmin, zmax = domain(distr)
	factor = 1
	r  = range(rmin, rmax, length=factor*(65))
	z  = range(zmin, zmax, length=factor*(129))
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

# ╔═╡ 1b03e78f-fe9a-4a97-8ec4-f732d65244d4
let
	rmin, rmax, zmin, zmax = domain(distr)
	factor = 1
	r  = range(rmin, rmax, length=ceil(Int, factor*65))
	z  = range(zmin, zmax, length=ceil(Int, factor*129))
	_I(r,z) = I(distr, r, z)
	nr, nz = length(r)-1, length(z)-1
	matrix = zeros(nr, nz)
	for i in 1:length(r)-1, j in 1:length(z)-1
		matrix[i,j] = torroidal_segment_yield(_I, r[i:i+1], z[j:j+1])[1] 
	end
	total = sum(matrix) # , matrix   moment 0 - okay
	rmids = collect(0.5(r[1:end-1] .+ r[2:end]))
	zmids = collect(0.5(z[1:end-1] .+ z[2:end]))
	# size(matrix), size(zmids), size(transpose(zmids))
	repeat(transpose(zmids), nr, 1)
	r1 = sum(repeat(rmids, 1, nz) .* matrix) / total	
	z1 = sum(matrix .* repeat(transpose(zmids), nr, 1)) / total	
	md"""## Moments of neutron emission distribution

	#  | Name    |   Values  
 	-: | :-----: | :---------:
	0  | total   |  $total   
	1  | r1, z1  |  $(r1), $(z1)
	"""
end

# ╔═╡ f328d540-9e59-47dd-8169-9b022be9e7f6


# ╔═╡ ee031b63-17e6-4047-9a75-0bfe2e825127
l

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
# ╠═a1070c6d-3bb2-4cd8-9c49-fc2a079814c2
# ╠═cbb99594-30c3-4ccf-9251-677fc6755750
# ╠═1b03e78f-fe9a-4a97-8ec4-f732d65244d4
# ╠═ad233063-b49b-49f2-a6d9-511324b30486
# ╠═f328d540-9e59-47dd-8169-9b022be9e7f6
# ╠═ee031b63-17e6-4047-9a75-0bfe2e825127
