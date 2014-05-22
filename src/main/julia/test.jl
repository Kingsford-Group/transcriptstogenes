require("argparse")
using ArgParse
#using Gadfly
using Gaston

function parseCommandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--input"
            help = "input file"
    end
    
    return parse_args(s)
end

function drawPlot(infile, outfile)
  Gaston.set_terminal("aqua")
  inputFile = open(infile, "r")

    pid = Float64[]
    for line in eachline(inputFile)
        toks = split(line)
        percentID = float(split(strip(toks[4],','),'=')[2])
        push!(pid, percentID)
    end
    e, counts = hist(pid)
    # draw(SVG(outfile, 6inch, 3inch), 
    #      plot(x=e, y=counts, Guide.XLabel("% identity"), Guide.YLabel("# of pairs"), Geom.normbar))
    histogram(pid, "bins", 100, "norm", 1, "color", "blue", "title", "redundancy")
    close(inputFile)
end

function main()
    parsedArgs = parseCommandline()
    println("Parsed args:")
    for pa in parsedArgs
        println("$(pa[1]) => $(pa[2])")
    end
    println(parsedArgs["input"])
    drawPlot(parsedArgs["intput"], "sample.svg")
        
end

#main()
