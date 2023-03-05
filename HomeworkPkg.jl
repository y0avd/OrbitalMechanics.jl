#=
Author: Yoav Dekel
=#

# homework submission functions
@enum baseUnit SEC KM KG RAD OTHER

# this function prints answers in red text, supports MathJax for Markdown in str input is raw"str"
function PrettyAns(str::String, ans::Union{String, Number}; unit::baseUnit=OTHER, unit_str::String="", sigfigs::Integer=3, ref::Bool=false)
    if typeof(ans) == String
        ans_str = ans
    else
        if unit == OTHER || length(unit_str) > 0
            ans_str = string(round(ans,digits=sigfigs))
        else
            ans_str = PrettyUnit(ans,unit,sigfigs=sigfigs)
        end
    end

    if ref
        display("text/markdown", str*" "*ans_str*unit_str)
    else
        display("text/markdown", str*" <b style=color:red;>"*ans_str*unit_str*"</b>")
    end
end

function PrettyUnit(ans::Number,unit::baseUnit;sigfigs::Integer=3)
    if unit == SEC
        yr = 365.24219
        if abs(ans)/60 < 1
            return string(round(ans,digits=sigfigs))*" seconds"
        elseif abs(ans)/(3600) < 1
            return string(round(ans/60,digits=sigfigs))*" minutes"
        elseif abs(ans)/(3600*24) < 1
            return string(round(ans/3600,digits=sigfigs))*" hours"
        elseif abs(ans)/(3600*24*yr) < 1 
            return string(round(ans/(24*3600),digits=sigfigs))*" days"
        else
            return string(round(ans/(24*3600*yr),digits=sigfigs))*" years"
        end
    elseif unit == KM
        AU = 1.495978707e8
        if abs(ans) < 1
            return string(round(1e3ans,digits=sigfigs))*" m"
        elseif abs(ans)/AU < 1
            return string(round(ans,digits=sigfigs))*" km"
        else
            return string(round(ans/AU,digits=sigfigs))*" AU"
        end
    elseif unit == KG
        if abs(ans) < 1
            return string(round(1e3ans,digits=sigfigs))*" g"
        else
            return string(round(ans,digits=sigfigs))*" kg"
        end
    elseif unit == RAD
        # always convert radians to degrees
        return string(round(rad2deg(ans),digits=sigfigs))*"°"
    end
end;