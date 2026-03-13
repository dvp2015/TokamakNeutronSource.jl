### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 26433bb2-6ff1-11ed-0943-01e79517b4f5
begin
	using Pkg
	Pkg.activate(temp=true)
	packages = [
		"Colors",
		"DataFrames",
		"GLMakie",
		"KernelDensity",
		"XLSX",
	]
	Pkg.add(packages)
	devdir = Pkg.devdir()
	eqdsk_pkg_dir = joinpath(devdir, "EQDSKReader.jl")
	tns_pkg_dir = joinpath(devdir, "TokamakNeutronSource.jl")
	Pkg.develop(path=eqdsk_pkg_dir)
	Pkg.develop(path=tns_pkg_dir)
	using EQDSKReader, TokamakNeutronSource
end

# ╔═╡ 9597cccc-5f54-4260-bf2e-2fcf4f3a169a
using 
	Colors,
	GLMakie,
    DataFrames,
	DelimitedFiles,
	KernelDensity,
	Printf,
	XLSX

# ╔═╡ cbf37834-2756-4138-9dc1-a679537c05a4
using	LinearAlgebra: norm


# ╔═╡ e8552241-d48d-4469-abae-484db725af36
using TokamakNeutronSource.PlasmaDistributions

# ╔═╡ 5f37ce53-477f-498d-b1ea-613fae26f3f8
begin
	eqdsk_path = joinpath(tns_pkg_dir, "test", "data", "beforeTQ.eqdsk")
	excel_path = joinpath(tns_pkg_dir, "test", "data", "TRT_215_8T_NBI.xlsx")
	@assert all(isfile, [eqdsk_path, excel_path])
end

# ╔═╡ 09221d19-80c5-42db-bc0c-e933bf0479f6
begin
	output_path = joinpath(tns_pkg_dir, "wrk", "plot-neutron-source-for-vant")
	mkpath(output_path)
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
	_I = PlasmaDistributions.I(distr, r, z)
	extrema(_I)
end

# ╔═╡ a2eda2d6-ddd8-4525-97b2-58112bcfb60b
md"""Total yield: $(@sprintf("%.2g", total_yield(distr)[1]))"""

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

# ╔═╡ e453e4a7-7fdb-4adb-9b31-d1f89f83bf13
md"""
Total yield and moment-0 of a distribution are the same.
"""

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

# ╔═╡ 4fd3535a-ed17-4c64-b8bb-a633c30e625f
begin
	open("../wrk/rbins.csv", "w") do io
		writedlm(io, rbins, ',')
	end
	open("../wrk/zbins.csv", "w") do io
		writedlm(io, zbins, ',')
	end
	open("../wrk/src.csv", "w") do io
		writedlm(io, src, ',')
	end
end

# ╔═╡ 6d363253-a395-4f6b-91e9-21e8f9fcabb8
maximum(src), sum(src), size(src), length(src)

# ╔═╡ b2df9bd8-64d6-4e86-9cbd-2b53e09cdc87
let
	ci = argmax(src)
	rbins[ci[1]], zbins[ci[2]]
end

# ╔═╡ 5d15e48f-4de8-4536-bf8a-8481862e6594
md"""
## Validate PTRAC data
"""

# ╔═╡ 67bb75ad-52ba-4e34-b7ce-00756c31e0de
datap = DelimitedFiles.readdlm("../wrk/ptrac.csv") .* 0.01  # cm -> m

# ╔═╡ adb6bdd8-e29d-4f80-8385-01e659b47387
let
	outliers = String[]
	foreach(datap[:,1], datap[:,2]) do r, z
		if 1.0 < ψ(distr)(r,z) 
			push!(outliers, "ψ($r, $z) = $(ψ(distr)(r,z))")
		end
	end
	if isempty(outliers)
		md"""All PTRAC points are within plasma volume."""
	else
		md"""
		$(length(outliers)) outliers:
		
		$(join(outliers, "\n\n"))
		"""
	end
end

# ╔═╡ a3a31b4b-2bd1-43e6-a4ad-97d99005b456
let
	f = Figure(size=(600, 800))
	ax = Axis(
		f[1,1]; 
		xlabel=L"R,m", 
		ylabel=L"Z,m", 
		aspect=DataAspect(),
		title=L"$$SDEF distribution (scatter)"
	)
	rmin, rmax, zmin, zmax = domain(distr)
	xlims!(ax, rmin, rmax)
	ylims!(ax, zmin, zmax)
	scatter!(
		ax, datap[:,1],
		datap[:,2],
		color=Colors.RGBA(0.25, 0.75, 0.35, 0.1),
		markersize=1
	)
	rbbs_points = [Point2f(x,y) for (x,y) in zip(eqdsk.rbbbs, eqdsk.zbbbs)]
	rlim_points = [Point2f(x,y) for (x,y) in zip(eqdsk.rlim, eqdsk.zlim)]
	scatter!(ax, eqdsk.rmaxis, eqdsk.zmaxis, color=:gray60, marker=:cross, label="Magnetic axis")
	lines!(ax, rbbs_points, label="Plasma boundary")
	lines!(ax, rlim_points, label="Limiter")
	axislegend(ax)
	# axs = Axis3(f[1,2])
	# surface!(axs, r,z,I)
	f	
end

# ╔═╡ b7d215b5-16c2-404b-8613-6d4814ca3beb
ptrac_kde = let
	rmin, rmax, zmin, zmax = domain(distr)
	np = 2^6
	kde(datap, npoints=(np, 2np), boundary=((rmin, rmax), (zmin, zmax)))
end

# ╔═╡ 62369e96-d0cf-4b61-8ee9-b6c049628179
let
	rmin, rmax, zmin, zmax = domain(distr)
	fontsize_theme = Theme(fontsize = 20, font = "JuliaMono")
	set_theme!(fontsize_theme, colormap = :bamako)
	f = Figure(;size=(700, 800))
	ax = Axis(
		f[1,1]; 
		xlabel=L"R,m", 
		ylabel=L"Z,m", 
		aspect=DataAspect(),
		title=L"$$SDEF distribution"
	)
	xlims!(ax, rmin, rmax)
	ylims!(ax, zmin, zmax)
	cntr = contourf!(ax,
		ptrac_kde,
		colormap=:bamako, 
		# colormap=:grays, 
		# linestyle="-",
		# linecolor=:black,
		# linewidth=2
	)
	rbbs_points = [Point2f(x,y) for (x,y) in zip(eqdsk.rbbbs, eqdsk.zbbbs)]
	rlim_points = [Point2f(x,y) for (x,y) in zip(eqdsk.rlim, eqdsk.zlim)]
	scatter!(ax, eqdsk.rmaxis, eqdsk.zmaxis, color=:gray60, marker=:cross, label="Magnetic axis")
	lines!(ax, rbbs_points, label="Plasma boundary")
	lines!(ax, rlim_points, label="Limiter")
	cb=Colorbar(
		f[1,2], 
		cntr, 
		label=L"$$Kernel Density", 
		# ticks=levels[2:2:end]
	)
	cb.tellheight = true
	cb.tellwidth = true
	# axislegend(ax)
	# axs = Axis3(f[1,3])
	# surface!(axs, r,z,ptrac_kde.density)
	f	
end

# ╔═╡ 94bfa226-cd13-4ea5-9ee3-017cd74c95d6
let
	rmin, rmax, zmin, zmax = domain(distr)
	x, y, z = ptrac_kde.x, ptrac_kde.y, ptrac_kde.density
	fontsize_theme = Theme(fontsize = 20, font = "JuliaMono")
	set_theme!(fontsize_theme, colormap = :bamako)
	f = Figure(size=(700, 800))
	ax = Axis3(
		f[1,1]; 
		xlabel=L"R,m", 
		ylabel=L"Z,m", 
		zlabel="SDEF pdf",
		aspect=(1, (zmax-zmin)/(rmax-rmin), 1),
		# title=L"$$SDEF distribution",
		azimuth=-3π/4,
		elevation = 1.2π/4,
		perspectiveness=0.5,
	)	
	# surf = surface!(ax, x, y, z, transparency=true, alpha=0.01)
	cntr = contour3d!(ax, x, y, z; 
		linewidth=3, levels = 20,
		# colormap = Reverse(:bamako) 
	)
	# wireframe!(ax, x, y, z, color=(:black, 0.5))

	# TODO dvp: the following causes bug
	# cb=Colorbar(
	# 	f[1,2], 
	# 	cntr, 
	# 	label=L"$$Kernel Density", 
	# 	# ticks=levels[2:2:end]
	# )
	f
end

# ╔═╡ e7b4f8b3-f3f0-4d5c-9080-eea93e01dd17
function plot_neutron_source(eqdsk::Content, distr::AbstractDistribution)
	f = Figure(; size=(600, 800))
	ax = Axis(
		f[1,1]; 
		xlabel=L"R,м", 
		ylabel=L"Z,м", 
		aspect=DataAspect(),
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
	)
	rbbs_points = [Point2f(x,y) for (x,y) in zip(eqdsk.rbbbs, eqdsk.zbbbs)]
	rlim_points = [Point2f(x,y) for (x,y) in zip(eqdsk.rlim, eqdsk.zlim)]
	scatter!(ax, eqdsk.rmaxis, eqdsk.zmaxis, color=:gray60, marker=:cross, label="Magnetic axis")
	lines!(ax, rbbs_points, label="Plasma boundary")
	lines!(ax, rlim_points, label="Limiter")
	cb=Colorbar(f[1,2], cntr, label=L"$I_{DD}$, $см^{-3}с^{-1}$", ticks=levels[1:2:end])
	cb.tellheight = true
	cb.tellwidth = true
	f	
end

# ╔═╡ 6888f022-b9dc-43f6-8759-772f97afc734
begin
	fig2 = plot_neutron_source(eqdsk, distr)
	save(joinpath(output_path, "fig2-dd-neutron-source.png"), fig2)
	fig2
end

# ╔═╡ 2b5b7b8f-6186-4e8b-b5a1-c5adac250f03
function ptrac_to_csv(dst, src)
	PREV_LINE = "       9000         1        40        47         0"
	is_data_line = false
	open(dst, "w") do out
		for line in readlines(src)
			if is_data_line
				x, y, z = map(x -> parse(Float64, x), split(line))
				r = LinearAlgebra.norm([x, y])
				println(out, r, ' ', z)
				is_data_line = false
			else
				is_data_line =  startswith(line, PREV_LINE)
			end
		end
	end
end


# ╔═╡ 13d29728-e0e3-4118-b73e-8509f798047d
isfile("../wrk/ptrac.csv") || ptrac_to_csv("../wrk/ptrac.csv", "../wrk/ptrac") 

# ╔═╡ a368bb91-e7ae-4fbc-85fe-ede9d95d01dc
md"""
## Plot I(ψ)
"""

# ╔═╡ 6c0c235d-8291-48ed-87c8-c0d64b161150
begin
	psi = excel[!, "ψ"]
	Y = I(distr)(psi)
end

# ╔═╡ 0f96bad4-4690-4401-b377-f1b792ba8add
insertcols!(excel, :n, (:I=>Y), after=true)

# ╔═╡ 52a270ab-ac9b-4dd0-b727-93123bc6d0ff
excel

# ╔═╡ 9e191d8e-4558-4cc2-8bf2-ae9ac70dd299
XLSX.writetable(joinpath(output_path, "psi-t-n-I.xlsx"), "ψ-distributions"=>excel, overwrite=true)

# ╔═╡ 350df4d4-d70e-4961-8ae0-8b8f602d89a3
function plot_neutron_source_vs_psi(excel::DataFrame)
	f = Figure()
	ax1 = Axis(
		f[1,1]; 
		ylabel=L"T_{i},\ keV",
		ylabelcolor=:blue,
		yticklabelcolor=:blue,
		limits=(0.0, 1.0, 0.0, 23),
		# xgridvisible=false,
		ygridvisible=false,
	)
	ax2 = Axis(
		f[1,1]; 
		ylabel=L"n_{i},\ 10^{14}cm^{-3}",
		ylabelcolor=:red,
		yticklabelcolor=:red,
		yaxisposition = :right,
		limits=(0.0, 1.0, 0.0, 1),
		xgridvisible=false,
		ygridvisible=false,
	)
	hidespines!(ax2)
	hidexdecorations!(ax2)
	ax3 = Axis(
		f[2,1]; 
		xlabel=L"ψ", 
		ylabel=L"I,\ 10^{10}cm^{-3}s^{-1}",
		ylabelcolor=:green,
		yticklabelcolor=:green,
		limits=(0.0, 1.0, 0, 1.5),
	)
	plot!(ax1, excel.:ψ, excel.:T, color=:blue )
	plot!(ax2, excel.:ψ, 1e-14*excel.:n, color=:red )
    plot!(ax3, excel.:ψ, 1e-10*excel.:I, color=:green )
	f	
end



# ╔═╡ d8e24d69-c661-47d6-9016-f07ca8272f36


# ╔═╡ 61954cdc-6ea6-45cd-979a-7add00c67208
fig3 = plot_neutron_source_vs_psi(excel)

# ╔═╡ 8625e77c-7941-4a99-9104-6f1ae5625ba7
save(joinpath(output_path, "fig3-neutron-source-vs-psi.png"), fig3)

# ╔═╡ Cell order:
# ╠═26433bb2-6ff1-11ed-0943-01e79517b4f5
# ╠═9597cccc-5f54-4260-bf2e-2fcf4f3a169a
# ╠═cbf37834-2756-4138-9dc1-a679537c05a4
# ╠═5f37ce53-477f-498d-b1ea-613fae26f3f8
# ╠═09221d19-80c5-42db-bc0c-e933bf0479f6
# ╠═50b9960f-7673-4975-b24b-d96ebb16d50f
# ╠═e8552241-d48d-4469-abae-484db725af36
# ╠═7e316d24-32ca-4084-b7ec-fa0b3cee1894
# ╠═11d8369b-a698-4cdb-96ed-5db4244f099e
# ╠═0f292ede-4a1f-4387-9f49-e0ac8ee993de
# ╟─a2eda2d6-ddd8-4525-97b2-58112bcfb60b
# ╠═6888f022-b9dc-43f6-8759-772f97afc734
# ╟─cbb99594-30c3-4ccf-9251-677fc6755750
# ╟─e453e4a7-7fdb-4adb-9b31-d1f89f83bf13
# ╟─ad233063-b49b-49f2-a6d9-511324b30486
# ╟─e2dde790-8e0f-4189-9971-a366ff2d6526
# ╠═e8381d26-3899-4f63-a905-aa7a47333401
# ╠═a952abea-f29e-46bf-9e7e-8e9bd16b1b2c
# ╟─6eec651c-e7ee-40e6-89f7-0b37af8e6c15
# ╠═187df34e-77a7-4c0f-9d10-ba1edfbdd41f
# ╠═1f5d7e83-7a13-4046-a7cc-7aa03e293d18
# ╠═216554fb-34d9-4ab8-83e3-bab6f438c5df
# ╠═4fd3535a-ed17-4c64-b8bb-a633c30e625f
# ╠═6d363253-a395-4f6b-91e9-21e8f9fcabb8
# ╠═b2df9bd8-64d6-4e86-9cbd-2b53e09cdc87
# ╟─5d15e48f-4de8-4536-bf8a-8481862e6594
# ╠═13d29728-e0e3-4118-b73e-8509f798047d
# ╟─67bb75ad-52ba-4e34-b7ce-00756c31e0de
# ╠═adb6bdd8-e29d-4f80-8385-01e659b47387
# ╠═a3a31b4b-2bd1-43e6-a4ad-97d99005b456
# ╠═b7d215b5-16c2-404b-8613-6d4814ca3beb
# ╠═62369e96-d0cf-4b61-8ee9-b6c049628179
# ╠═94bfa226-cd13-4ea5-9ee3-017cd74c95d6
# ╠═e7b4f8b3-f3f0-4d5c-9080-eea93e01dd17
# ╠═2b5b7b8f-6186-4e8b-b5a1-c5adac250f03
# ╠═a368bb91-e7ae-4fbc-85fe-ede9d95d01dc
# ╠═6c0c235d-8291-48ed-87c8-c0d64b161150
# ╠═0f96bad4-4690-4401-b377-f1b792ba8add
# ╠═52a270ab-ac9b-4dd0-b727-93123bc6d0ff
# ╠═9e191d8e-4558-4cc2-8bf2-ae9ac70dd299
# ╠═350df4d4-d70e-4961-8ae0-8b8f602d89a3
# ╠═d8e24d69-c661-47d6-9016-f07ca8272f36
# ╠═61954cdc-6ea6-45cd-979a-7add00c67208
# ╠═8625e77c-7941-4a99-9104-6f1ae5625ba7
