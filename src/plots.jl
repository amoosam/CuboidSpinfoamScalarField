ENV["GKSwstype"] = "nul"
using Plots    #uses GR Backend

function plot_acc(alpha::Float64, mass::Float64, acc::Array{<:AbstractFloat,1},filenum)
    therms = 1:size(acc)[1]
    plot(therms,acc./10000, seriestype= :scatter, title="alpha: $alpha, mass: $mass", titleloc= :right, titlefontsize=14)
    xlabel!("Thermalization steps (x10,000)")
    ylabel!("acceptance rate per 10,000")
    savefig(string(path,"therm_acc-",filenum,ver,".png"))
end

function plot_vol(vol::Array{<:AbstractFloat,1},filenum)
    therms = 1:size(vol)[1]
    plot(therms[10000:end],log.(vol[10000:end]))
    xlabel!("Thermalization steps")
    ylabel!("log(4-volume)")
    savefig(string(path,"therm_vol-",filenum,ver,".png"))
end
function plot_amp(ampl::Array{<:AbstractFloat,1},filenum)
    therms = 1:size(ampl)[1]
    plot(therms[10000:end],log.(10,ampl[10000:end]))
    xlabel!("thermalization steps")
    ylabel!("log(amplitude)")
    savefig(string(path,"therm_ampl-",filenum,ver,".png"))
end

function plot_thermresults(alpha::Float64, mass::Float64, ampl::Array{<:AbstractFloat,1},vol::Array{<:AbstractFloat,1}, filenum::String)
    l = @layout[a;b]
    therms = 1:size(vol)[1]
    #p1 = plot(therms[10000:end],log.(10,vol[10000:end]), yscale = :log10, title="alpha: $alpha, mass: $mass", titleloc= :right, titlefontsize=12, legend=false)
    p1 = plot(therms[10000:end],vol[10000:end], yscale = :log10, title="alpha: $alpha, mass: $mass", titleloc= :right, titlefontsize=12, legend=false)
    xlabel!("Thermalization steps")
    ylabel!("log(4-volume)")

    #p2 = plot(therms[10000:end],log.(10,ampl[10000:end]), yscale=:log10, legend=false)
    p2 = plot(therms[10000:end],ampl[10000:end], yscale=:log10, legend=false)
    xlabel!("thermalization steps")
    ylabel!("log(amplitude)")

    plot(p1,p2, layout=l)
    savefig(string(path,"therm_results-",filenum,ver,".png"))
end

function plot_therm_acc_results(alpha::Float64, mass::Float64, acc::Array{<:AbstractFloat,1},ampl::Array{<:AbstractFloat,1},vol::Array{<:AbstractFloat,1}, filename::String)
    l = @layout[a;b;c]

    therms = 1:size(vol)[1]
    p1 = plot(therms[2:end],log.(10,vol[2:end]), title="alpha: $alpha, mass: $mass", titleloc= :right, titlefontsize=12, legend=false)
    xlabel!("Thermalization steps")
    ylabel!("log(4-volume)")

    #p2 = plot(therms[2:end],log.(10,ampl[2:end]), yscale = :log10, legend=false)
    p2 = plot(therms[2:end],log.(10,ampl[2:end]), legend=false)
    xlabel!("thermalization steps")
    ylabel!("log(amplitude)")

    therms = 1:size(acc)[1]
    p3 = plot(therms,acc, seriestype= :scatter, title="alpha: $alpha, mass:$mass", titleloc= :right, titlefontsize=14)
    xlabel!("Thermalization steps (x10,000)")
    ylabel!("acceptance rate per 10,000")

    plot(p1,p2,p3, layout=l, size=(900,600))
    savefig(string(path,"therm_acc_results-",filename,ver,".png"))
end

function plot_results(alpha::Float64, mass::Float64, acc::Array{<:AbstractFloat,1},ampl::Array{<:AbstractFloat,1}, filename::String)
    l = @layout[a;b]
    
    samples = 1:size(acc)[1]
    p1 = plot(samples,acc, title="alpha: $alpha, mass: $mass", titleloc= :right, titlefontsize=12, legend=false )
    xlabel!("samples")
    ylabel!("acceptance rate")

    p2 = plot(samples,log.(10,ampl), legend=false)
    xlabel!("samples")
    ylabel!("log(amplitude)")

    plot(p1,p2, layout=l)
    savefig(string(path,"samp_results-",filename,ver,".png"))
end


