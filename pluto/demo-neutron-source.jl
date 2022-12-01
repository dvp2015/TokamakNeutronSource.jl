### A Pluto.jl notebook ###
# v0.19.16

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
using TokamakNeutronSource.PlasmaDistribution

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

# ╔═╡ f453daec-8564-48e4-8869-834d3ef64fd7
let
	_ψ = distr.ψ
	_I = I(distr)
	__I(r,z) = _I(_ψ(r,z))
end

# ╔═╡ 2c9df88e-b0a1-4078-a539-685c3c85f676
function compute_total_output(distr)
	R0 = distr.rmin
	ΔR = distr.rmax - R0
	Z0 = distr.zmin
	ΔZ = distr.zmax - Z0
	_ψ = distr.ψ
	_I = I(distr)
	__I(r,z) = _I(_ψ(r,z))
	function integrand(x, f) 
		r = ΔR * x[1] + R0
		z = ΔZ * x[2] + Z0
		f[1] = r*__I(r, z)
	end
	# cuhre(integrand, rtol=1e-4)
	integral, error, probability, neval, fail, nregions = cuhre(integrand)
	# Normalization:
	# 2π - integral over torroidal direction, 
	# (...)(...) - area of R,Z integration domain, square meters
	# 1e6 - m^3/cm^3
	normalization = 2.0e6π*ΔR*ΔZ
	(
		integral[1]*normalization,
		error[1]*normalization,
		neval,
		fail
	)
end

# ╔═╡ a3e884fe-362a-4be6-aba6-6cd33da663ce
@time compute_total_output(distr)

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
# ╠═f453daec-8564-48e4-8869-834d3ef64fd7
# ╠═2c9df88e-b0a1-4078-a539-685c3c85f676
# ╠═a3e884fe-362a-4be6-aba6-6cd33da663ce
