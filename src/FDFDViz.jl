module FDFDViz

using FDFD, Plots, AxisArrays, Printf

export plot_field, plot_device

"    plot_field(ax, field::Field; cbar::Bool=false, funcz=real)"
function plot_field(field::Field; cbar::Bool=false, funcz=real, component=nothing)
	# Set defaults
	isa(field, FieldTM) && component === nothing && (component=:Ez);
	isa(field, FieldTE) && component === nothing && (component=:Hz);
	isa(field, FieldAll) && component === nothing && (component=:Ex);
	
	Z = funcz.(permutedims(field[:, :, component], (2,1)))

	if funcz == abs
		vmin = 0;
		vmax = +maximum(abs.(Z));
		cmap = :magma;
	elseif funcz == real
		Zmx  = maximum(abs.(Z));
		vmin = -maximum(Zmx);
		vmax = +maximum(Zmx);
		cmap = :RdBu;
	else
		error("Unknown function specified.");
	end

	extents = [ field.grid.bounds[1][1], field.grid.bounds[2][1],
	            field.grid.bounds[1][2], field.grid.bounds[2][2] ];

	p = heatmap(Z, c=cmap, extent=extents, origin="lower", clims=(vmin, vmax), label="||E||")
	xlabel!("x");
	ylabel!("y");
	title!("ω/2π = $(real(field.ω/2π/1e12)) THz");

	return p
end

"    plot_device(device::AbstractDevice; outline::Bool=false, lc::String=\"k\", lcm::String=\"k\")"
function plot_device(device::AbstractDevice; outline::Bool=false, lc::String="k", lcm::String="k")
	return plot_device(device; outline=outline, lc=lc, lcm=lcm);
end

"    plot_device(ax, device::AbstractDevice; outline::Bool=true, lc::String=\"k\", lcm::String=\"k\")"
function plot_device(ax, device::AbstractDevice; outline::Bool=true, lc::String="k", lcm::String="k")
	Z = real.(permutedims(device.ϵᵣ, (2,1)));
	Zi = imag.(permutedims(device.ϵᵣ, (2,1)));

	p = plot()

	if outline
		contour!(xc(device.grid), yc(device.grid), Z, linewidths=0.25, colors=lc);
	else
		extents = [ device.grid.bounds[1][1], device.grid.bounds[2][1],
		            device.grid.bounds[1][2], device.grid.bounds[2][2] ];
		heatmap!(Z, extent=extents, origin="lower", c=:YlGnBu);
	end

	contour!(xc(device.grid), yc(device.grid), Zi, colors=lc, linewidths=0.5, linestyles=":");

	if isa(device, ModulatedDevice)
		Z2 = abs.(permutedims(device.Δϵᵣ, (2,1)));
		contour!(xc(device.grid), yc(device.grid), Z2, levels=1, linewidths=0.5, colors=lcm);
	end
	xlabel!("x");
	ylabel!("y");

	return p
end

end
