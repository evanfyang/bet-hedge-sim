import JSON
export read_parameters

# read parameters from JSON file
function read_parameters(defpars, inputfile=nothing)
    pars = copy(defpars)

    # read JSON file
    if inputfile != nothing
        inpars = JSON.parsefile(inputfile)
    else
        inpars = Dict()
    end

    for parkey in keys(defpars)
        if "type" in keys(pars[parkey]) &&
            isprimitivetype(pars[parkey]["type"])
            T = pars[parkey]["type"]
        else
            # Default type is Float
            T = Float64
        end
        if T <: Int
            convertf = (val) -> round(T, val)
        else
            convertf = (val) -> convert(T, val)
        end

        # use defpars for list of usabe parameters in JSON
        if parkey in keys(inpars)
            if "value" in keys(inpars[parkey])
                val = inpars[parkey]["value"]
            elseif "range" in keys(inpars[parkey])
                valr = inpars[parkey]["range"]
                if "log" in keys(valr)
                    b = valr["log"]
                    rf = (r) -> b .^ r
                    pop!(valr, "log")
                else
                    rf = (r) -> r
                end
                start = pop!(valr, "start")
                rkws = Dict(zip(Symbol.(keys(valr)), values(valr)))
                val = rf(range(start; rkws...))
            end
        else
            val = pars[parkey]["value"]
        end

        if !isstructtype(typeof(val))
            pars[parkey] = [convertf(val)]
        else
            pars[parkey] = convertf.(val)
        end
    end

    return pars
end